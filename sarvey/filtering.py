#!/usr/bin/env python

# SARvey - A multitemporal InSAR time series tool for the derivation of displacements.
#
# Copyright (C) 2021-2024 Andreas Piter (IPI Hannover, piter@ipi.uni-hannover.de)
#
# This software was developed together with FERN.Lab (fernlab@gfz-potsdam.de) in the context
# of the SAR4Infra project with funds of the German Federal Ministry for Digital and
# Transport and contributions from Landesamt fuer Vermessung und Geoinformation
# Schleswig-Holstein and Landesbetrieb Strassenbau und Verkehr Schleswig-Holstein.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# Important: This package uses PyMaxFlow. The core of PyMaxflows library is the C++
# implementation by Vladimir Kolmogorov. It is also licensed under the GPL, but it REQUIRES that you
# cite [BOYKOV04] (see LICENSE) in any resulting publication if you use this code for research purposes.
# This requirement extends to SARvey.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.

"""Filtering module for SARvey."""
import time
import multiprocessing
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
import gstools as gs
from logging import Logger

from mintpy.utils import ptime

import sarvey.utils as ut


def launchSpatialFiltering(parameters: tuple):
    """Launch_spatial_filtering.

    Launches the spatial filtering to estimate the atmospheric phase screen with low-pass filtering.

    Parameters
    ----------
    parameters: tuple
        Tuple containing the following parameters:

        idx_range: np.ndarray
            range of indices for the time series
        num_time: int
            number of time steps
        residuals: np.ndarray
            residual phase (size: num_points x num_ifgs)
        coord_utm1: np.ndarray
            coordinates in UTM of the first-order points for which the residuals are given (size: num_points_p1 x 2)
        coord_utm2: np.ndarray
            coordinates in UTM of the new points which shall be interpolated (size: num_points_p2 x 2)
        bins: np.ndarray
            bin edges for the variogram
        bool_plot: bool
            boolean flag to plot intermediate results
        logger: Logger
            Logging handler

    Returns
    -------
    idx_range: np.ndarray
        range of indices for the time series
    aps1: np.ndarray
        atmospheric phase screen for the known points (size: num_points_p1 x num_ifgs)
    aps2: np.ndarray
        atmospheric phase screen for the new points (size: num_points_p2 x num_ifgs)
    """
    # Unpack the parameters
    (idx_range, num_time, residuals, coord_utm1, coord_utm2, bins, bool_plot, logger) = parameters

    x = coord_utm1[:, 1]
    y = coord_utm1[:, 0]
    x_new = coord_utm2[:, 1]
    y_new = coord_utm2[:, 0]

    aps1 = np.zeros((coord_utm1.shape[0], num_time), dtype=np.float32)
    aps2 = np.zeros((coord_utm2.shape[0], num_time), dtype=np.float32)

    prog_bar = ptime.progressBar(maxValue=num_time)

    for i in range(num_time):
        field = residuals[:, i].astype(np.float32)

        # 1) estimate the variogram of the field
        bin_center, vario = gs.vario_estimate(pos=[x, y], field=field, bin_edges=bins)

        # 2) fit model to empirical variogram
        model = gs.Stable(dim=2)
        try:
            model.fit_variogram(x_data=bin_center, y_data=vario, nugget=True, max_eval=1500)
        except RuntimeError as err:
            logger.error(msg="\nIMAGE {}: Not able to fit variogram! {}".format(idx_range[i], err))
            if bool_plot:
                fig, ax = plt.subplots(2, figsize=[10, 5])
                sca1 = ax[0].scatter(x, y, c=field)
                plt.colorbar(sca1, ax=ax[0], pad=0.03, shrink=0.5)
                ax[0].set_title("Not able to fit variogram! - PS1 residuals")
                ax[1].scatter(bin_center, vario)
                ax[1].set_xlabel("distance in [m]")
                ax[1].set_ylabel("semi-variogram")
                plt.close(fig)
            prog_bar.update(value=i + 1, every=1, suffix='{}/{} images'.format(i + 1, num_time))
            continue

        # 3) estimate parameters of kriging
        sk = gs.krige.Simple(
            model=model,
            cond_pos=[x, y],
            cond_val=field,
        )

        # 4) evaluate the kriging model at ORIGINAL locations
        fld_sk, _ = sk((x, y), return_var=True)
        aps1[:, i] = fld_sk

        # 5) evaluate the kriging model at NEW locations
        fld_sk_new, var_sk_new = sk((x_new, y_new), return_var=True)
        aps2[:, i] = fld_sk_new

        prog_bar.update(value=i + 1, every=1, suffix='{}/{} images'.format(i + 1, num_time))

        # 5) show results
        if bool_plot:
            min_val = np.min(field)
            max_val = np.max(field)

            fig, ax = plt.subplots(2, 2, figsize=[10, 5])

            cur_ax = ax[0, 0]
            sca1 = cur_ax.scatter(x, y, c=field, vmin=min_val, vmax=max_val)
            plt.colorbar(sca1, ax=cur_ax, pad=0.03, shrink=0.5)
            cur_ax.set_title("PS1 residuals")

            cur_ax = ax[0, 1]
            cur_ax = model.plot(x_max=bin_center[-1], ax=cur_ax)
            cur_ax.scatter(bin_center, vario)
            cur_ax.set_xlabel("distance in [m]")
            cur_ax.set_ylabel("semi-variogram")

            if coord_utm2 is not None:
                cur_ax = ax[1, 0]
                sca2 = cur_ax.scatter(x_new, y_new, c=fld_sk_new, vmin=min_val, vmax=max_val)
                plt.colorbar(sca2, ax=cur_ax, pad=0.03, shrink=0.5)
                cur_ax.set_title("PS2 prediction of atmospheric effect")

                cur_ax = ax[0, 1]
                sca4 = cur_ax.scatter(x_new, y_new, c=var_sk_new)
                plt.colorbar(sca4, ax=cur_ax, pad=0.03, shrink=0.5)
                cur_ax.set_title("Variance of predicted atmospheric effect")

            plt.close(fig)

    return idx_range, aps1, aps2


def estimateAtmosphericPhaseScreen(*, residuals: np.ndarray, coord_utm1: np.ndarray, coord_utm2: np.ndarray,
                                   num_cores: int = 1, bool_plot: bool = False,
                                   logger: Logger) -> tuple[np.ndarray, np.ndarray]:
    """Estimate_atmospheric_phase_screen.

    Estimates the atmospheric phase screen from a stack of phase time series for a sparse set of points.
    Kriging is used to estimate the spatial dependence and to interpolate the phase screen over a set of new points.

    Parameters
    ----------
    residuals: np.ndarray
        residual phase (size: num_points1 x num_images)
    coord_utm1: np.ndarray
        coordinates in UTM of the points for which the residuals are given (size: num_points1 x 2)
    coord_utm2: np.ndarray
        coordinates in UTM of the new points which shall be interpolated (size: num_points2 x 2)
    num_cores: int
        Number of cores
    bool_plot: bool
        boolean flag to plot intermediate results (default: False)
    logger: Logger
        Logging handler

    Returns
    -------
    aps1: np.ndarray
        atmospheric phase screen for the known points (size: num_points1 x num_images)
    aps2: np.ndarray
        atmospheric phase screen for the new points (size: num_points2 x num_images)
    """
    msg = "#" * 10
    msg += " ESTIMATE ATMOSPHERIC PHASE SCREEN (KRIGING) "
    msg += "#" * 10
    logger.info(msg=msg)

    start_time = time.time()

    num_points1 = residuals.shape[0]
    num_points2 = coord_utm2.shape[0]
    num_time = residuals.shape[1]  # can be either num_ifgs or num_images

    bins = gs.variogram.standard_bins(pos=(coord_utm1[:, 1], coord_utm1[:, 0]),
                                      dim=2, latlon=False, mesh_type='unstructured', bin_no=30, max_dist=None)

    if num_cores == 1:
        args = (np.arange(0, num_time), num_time, residuals, coord_utm1, coord_utm2, bins, bool_plot, logger)
        _, aps1, aps2 = launchSpatialFiltering(parameters=args)
    else:
        logger.info(msg="start parallel processing with {} cores.".format(num_cores))
        pool = multiprocessing.Pool(processes=num_cores)

        aps1 = np.zeros((num_points1, num_time), dtype=np.float32)
        aps2 = np.zeros((num_points2, num_time), dtype=np.float32)

        num_cores = num_time if num_cores > num_time else num_cores  # avoids having more samples than cores
        idx = ut.splitDatasetForParallelProcessing(num_samples=num_time, num_cores=num_cores)

        args = [(
            idx_range,
            idx_range.shape[0],
            residuals[:, idx_range],
            coord_utm1,
            coord_utm2,
            bins,
            False,
            logger) for idx_range in idx]

        results = pool.map(func=launchSpatialFiltering, iterable=args)

        # retrieve results
        for i, aps1_i, aps2_i in results:
            aps1[:, i] = aps1_i
            aps2[:, i] = aps2_i

    m, s = divmod(time.time() - start_time, 60)
    logger.debug(msg='time used: {:02.0f} mins {:02.1f} secs.\n'.format(m, s))

    return aps1, aps2


def simpleInterpolation(*, residuals: np.ndarray, coord_utm1: np.ndarray, coord_utm2: np.ndarray,
                        interp_method: str = "linear"):
    """SimpleInterpolation.

    Simple interpolation of atmospheric phase screen using scipy's griddata function with options "linear" or "cubic".
    For pixels outside the convex hull of the input points, the nearest neighbor is used.

    Parameters
    ----------
    residuals: np.ndarray
        residual phase (size: num_points x num_ifgs)
    coord_utm1: np.ndarray
        coordinates in UTM of the points for which the residuals are given (size: num_points_p1 x 2)
    coord_utm2: np.ndarray
        coordinates in UTM of the new points which shall be interpolated (size: num_points_p2 x 2)
    interp_method: str
        interpolation method (default: "linear"; options: "linear", "cubic")

    Returns
    -------
    aps1: np.ndarray
        atmospheric phase screen for the known points (size: num_points_p1 x num_images)
    aps2: np.ndarray
        atmospheric phase screen for the new points (size: num_points_p2 x num_images)
    """
    num_points2 = coord_utm2.shape[0]
    num_images = residuals.shape[1]

    aps1 = np.zeros_like(residuals, dtype=np.float32)
    aps2 = np.zeros((num_points2, num_images), dtype=np.float32)
    for i in range(num_images):
        aps1[:, i] = griddata(coord_utm1, residuals[:, i], coord_utm1, method=interp_method)
        aps2[:, i] = griddata(coord_utm1, residuals[:, i], coord_utm2, method=interp_method)
        # interpolation with 'linear' or 'cubic' yields nan values for pixel that need to be extrapolated.
        # interpolation with 'knn' solves this problem.
        mask_extrapolate = np.isnan(aps2[:, i])
        aps2[mask_extrapolate, i] = griddata(
            coord_utm1,
            residuals[:, i],
            coord_utm2[mask_extrapolate, :],
            method='nearest'
        )

    return aps1, aps2


def extractModelParams(model: gs.covmodel.models, logger: Logger):
    """
    Extract parameters from the given gstools model.

    Parameters
    ----------
    model : gs.covmodel.models
        The model from which to extract parameters.
    logger : Logger
        Logging handler.

    Returns
    -------
    tuple
        A tuple containing the model parameters.

    Raises
    ------
    ValueError
        If the model is not implemented.
    """
    if model.name == 'Stable':
        params = (model.var, model.len_scale, model.nugget, model.alpha)
        logger.debug(msg=f"Extract {len(params)} parameters from gs model {model.name}.")
    # elif model.name == 'Gaussian':
    #     params = (model.var, model.len_scale, model.nugget)
    else:
        error_msg = f"Model {model.name} not implemented yet."
        logger.error(msg=error_msg)
        raise ValueError(error_msg)
    return params


def applyModelParams(model: gs.covmodel.models, params: tuple, logger: Logger):
    """
    Apply parameters to the given gstools model.

    Parameters
    ----------
    model : gs.covmodel.models
        The model to which parameters will be applied.
    params : tuple
        A tuple containing the model parameters.
    logger : Logger
        Logging handler.

    Returns
    -------
    model
        The model with applied parameters.

    Raises
    ------
    ValueError
        If the model is not implemented.
    """

    if model.name == 'Stable':
        logger.debug(msg=f"Applying {params.size} parameters to gs model {model.name}.")
        model.var, model.len_scale, model.nugget, model.alpha = params
    # elif model.name == 'Gaussian':  # TODO: implement other models
    #     model.var, model.len_scale, model.nugget = params
    else:
        error_msg = f"Model {model.name} not implemented yet."
        logger.error(msg=error_msg)
        raise ValueError(error_msg)
    return model
