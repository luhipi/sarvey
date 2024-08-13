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

"""Unwrapping module for SARvey."""
import multiprocessing
from os.path import join, dirname
import time
from typing import Union

import matplotlib.pyplot as plt
import numpy as np
from kamui import unwrap_arbitrary
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import structural_rank
from scipy.sparse.linalg import lsqr
from scipy.optimize import minimize
from logging import Logger

from mintpy.utils import ptime

import sarvey.utils as ut
from sarvey.ifg_network import IfgNetwork
from sarvey.objects import Network, NetworkParameter, AmplitudeImage


def objFuncTemporalCoherence(x, *args):
    """Compute temporal coherence from parameters and phase. To be used as objective function for optimization.

    Parameters
    ----------
    x: np.ndarray
        Search space for the DEM error in a 1D grid.
    args: tuple
        Additional arguments: (design_mat, obs_phase, scale_vel, scale_demerr).

    Returns
    -------
    1 - gamma: float
    """
    (design_mat, obs_phase, scale_vel, scale_demerr) = args

    # equalize the gradients in both directions
    x[0] *= scale_demerr
    x[1] *= scale_vel

    pred_phase = np.matmul(design_mat, x)
    res = (obs_phase - pred_phase.T).ravel()
    gamma = np.abs(np.mean(np.exp(1j * res)))
    return 1 - gamma


def gridSearchTemporalCoherence(*, demerr_grid: np.ndarray, vel_grid: np.ndarray, design_mat: np.ndarray,
                                obs_phase: np.ndarray):
    """Grid search which maximizes the temporal coherence as the objective function.

    Parameters
    ----------
    demerr_grid: np.ndarray
        Search space for the DEM error in a 2D grid.
    vel_grid: np.ndarray
        Search space for the velocity in a 2D grid.
    design_mat: np.ndarray
        Design matrix for estimating parameters from arc phase.
    obs_phase: np.ndarray
        Observed phase of the arc.

    Returns
    -------
    demerr: float
        estimated DEM error.
    vel: float
        estimated velocity.
    gamma: float
        estimated temporal coherence.
    """
    demerr_grid_flat = demerr_grid.flatten()
    vel_grid_flat = vel_grid.flatten()
    gamma_flat = np.array(
        [1 - objFuncTemporalCoherence(np.array([demerr_grid_flat[i], vel_grid_flat[i]]),
                                      design_mat, obs_phase, 1, 1)
         for i in range(demerr_grid_flat.shape[0])])
    gamma = gamma_flat.reshape(demerr_grid.shape)
    idx_max_gamma = np.argmax(gamma_flat)

    # return demerr_grid_flat[idx_max_gamma], vel_grid_flat[idx_max_gamma], gamma_flat[idx_max_gamma]
    return demerr_grid_flat[idx_max_gamma], vel_grid_flat[idx_max_gamma], gamma


def findOptimum(*, obs_phase: np.ndarray, design_mat: np.ndarray, val_range: np.ndarray):
    """Find optimal value within a one dimensional search space that fits to the observed phase.

    Parameters
    ----------
    obs_phase: np.ndarray
        Observed phase of the arc.
    design_mat: np.ndarray
        Design matrix for estimating parameters from arc phase.
    val_range: np.ndarray
        Range of possible values for the solution. Can be either for DEM error or velocity.

    Returns
    -------
    opt_val: scipy.optimize.minimize return value
    gamma: float
    pred_phase: np.ndarray
    """
    pred_phase = design_mat[:, np.newaxis] * val_range[np.newaxis, :]  # broadcasting
    if len(obs_phase.shape) == 2:
        # step densification
        res = obs_phase[:, np.newaxis, :] - pred_phase.T
        res = np.moveaxis(res, 0, 1)
        res = res.reshape((pred_phase.shape[1], -1))  # combine residuals from all arcs
    else:
        # step consistency check
        res = obs_phase - pred_phase.T

    gamma = np.abs(np.mean(np.exp(1j * res), axis=1))
    max_idx = np.argmax(gamma)
    opt_val = val_range[max_idx]
    return opt_val, gamma[max_idx], pred_phase[:, max_idx]


def oneDimSearchTemporalCoherence(*, demerr_range: np.ndarray, vel_range: np.ndarray, obs_phase: np.ndarray,
                                  design_mat: np.ndarray):
    """One dimensional search for maximum temporal coherence that fits the observed arc phase.

    Parameters
    ----------
    demerr_range: np.ndarray
        Search space for the DEM error in a 1D grid.
    vel_range: np.ndarray
        Search space for the velocity in a 1D grid.
    design_mat: np.ndarray
        Design matrix for estimating parameters from arc phase.
    obs_phase: np.ndarray
        Observed phase of the arc.

    Returns
    -------
    demerr: float
    vel: float
    gamma: float
    """
    demerr, gamma_demerr, pred_phase_demerr = findOptimum(
        obs_phase=obs_phase,
        design_mat=design_mat[:, 0],
        val_range=demerr_range
    )

    vel, gamma_vel, pred_phase_vel = findOptimum(
        obs_phase=obs_phase,
        design_mat=design_mat[:, 1],
        val_range=vel_range
    )

    if gamma_vel > gamma_demerr:
        demerr, gamma_demerr, pred_phase_demerr = findOptimum(
            obs_phase=obs_phase - pred_phase_vel,
            design_mat=design_mat[:, 0],
            val_range=demerr_range
        )
        vel, gamma_vel, pred_phase_vel = findOptimum(
            obs_phase=obs_phase - pred_phase_demerr,
            design_mat=design_mat[:, 1],
            val_range=vel_range
        )
    else:
        vel, gamma_vel, pred_phase_vel = findOptimum(
            obs_phase=obs_phase - pred_phase_demerr,
            design_mat=design_mat[:, 1],
            val_range=vel_range
        )
        demerr, gamma_demerr, pred_phase_demerr = findOptimum(
            obs_phase=obs_phase - pred_phase_vel,
            design_mat=design_mat[:, 0],
            val_range=demerr_range
        )

    # improve initial estimate with gradient descent approach
    scale_demerr = demerr_range.max()
    scale_vel = vel_range.max()

    demerr, vel, gamma = gradientSearchTemporalCoherence(
        scale_vel=scale_vel,
        scale_demerr=scale_demerr,
        obs_phase=obs_phase,
        design_mat=design_mat,
        x0=np.array([demerr / scale_demerr,
                     vel / scale_vel]).T
    )

    pred_phase = np.matmul(design_mat, np.array([demerr, vel]))
    res = (obs_phase - pred_phase.T).ravel()
    gamma = np.abs(np.mean(np.exp(1j * res)))
    return demerr, vel, gamma


def gradientSearchTemporalCoherence(*, scale_vel: float, scale_demerr: float, obs_phase: np.ndarray,
                                    design_mat: np.ndarray, x0: np.ndarray):
    """GradientSearchTemporalCoherence.

    Parameters
    ----------
    scale_demerr: float
        Scaling factor for DEM error to equalize the axis of the search space.
    scale_vel: float
        Scaling factor for velocity to equalize the axis of the search space.
    design_mat: np.ndarray
        Design matrix for estimating parameters from arc phase.
    obs_phase: np.ndarray
        Observed phase of the arc.
    x0: np.ndarray
        Initial values for optimization.

    Returns
    -------
    demerr: float
    vel: float
    gamma: float
    """
    opt_res = minimize(
        objFuncTemporalCoherence,
        x0,
        args=(design_mat, obs_phase, scale_vel, scale_demerr),
        bounds=((-1, 1), (-1, 1)),
        method='L-BFGS-B'
    )
    gamma = 1 - opt_res.fun
    demerr = opt_res.x[0] * scale_demerr
    vel = opt_res.x[1] * scale_vel
    return demerr, vel, gamma


def launchAmbiguityFunctionSearch(parameters: tuple):
    """Wrap for launching ambiguity function for temporal unwrapping in parallel.

    Parameters
    ----------
    parameters: tuple
        Arguments for temporal unwrapping in parallel.

    Returns
    -------
    arc_idx_range: np.ndarray
    demerr: np.ndarray
    vel: np.ndarray
    gamma: np.ndarray
    """
    (arc_idx_range, num_arcs, phase, slant_range, loc_inc, ifg_net_obj, wavelength, velocity_bound, demerr_bound,
     num_samples) = parameters

    demerr = np.zeros((num_arcs, 1), dtype=np.float32)
    vel = np.zeros((num_arcs, 1), dtype=np.float32)
    gamma = np.zeros((num_arcs, 1), dtype=np.float32)

    design_mat = np.zeros((ifg_net_obj.num_ifgs, 2), dtype=np.float32)

    demerr_range = np.linspace(-demerr_bound, demerr_bound, num_samples)
    vel_range = np.linspace(-velocity_bound, velocity_bound, num_samples)

    # prog_bar = ptime.progressBar(maxValue=num_arcs)

    factor = 4 * np.pi / wavelength

    for k in range(num_arcs):
        design_mat[:, 0] = factor * ifg_net_obj.pbase_ifg / (slant_range[k] * np.sin(loc_inc[k]))
        design_mat[:, 1] = factor * ifg_net_obj.tbase_ifg

        demerr[k], vel[k], gamma[k] = oneDimSearchTemporalCoherence(
            demerr_range=demerr_range,
            vel_range=vel_range,
            obs_phase=phase[k, :],
            design_mat=design_mat
        )

    return arc_idx_range, demerr, vel, gamma


def temporalUnwrapping(*, ifg_net_obj: IfgNetwork, net_obj: Network,  wavelength: float, velocity_bound: float,
                       demerr_bound: float, num_samples: int, num_cores: int = 1, logger: Logger) -> \
        tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Solve ambiguities for every arc in spatial Network object.

    Parameters
    ----------
    ifg_net_obj: IfgNetwork
        The IfgNetwork object.
    net_obj: Network
        The Network object.
    wavelength: float
        The wavelength.
    velocity_bound: float
        The velocity bound.
    demerr_bound: float
        The DEM error bound.
    num_samples: int
        The number of samples for the search space.
    num_cores: int
        Number of cores to be used. Default is 1.
    logger: Logger
        Logging handler.

    Returns
    -------
    demerr: np.ndarray
    vel: np.ndarray
    gamma: np.ndarray
    """
    msg = "#" * 10
    msg += " TEMPORAL UNWRAPPING: AMBIGUITY FUNCTION "
    msg += "#" * 10
    logger.info(msg)

    start_time = time.time()

    if num_cores == 1:
        args = (
            np.arange(net_obj.num_arcs), net_obj.num_arcs, net_obj.phase,
            net_obj.slant_range, net_obj.loc_inc, ifg_net_obj, wavelength, velocity_bound, demerr_bound, num_samples)
        arc_idx_range, demerr, vel, gamma = launchAmbiguityFunctionSearch(parameters=args)
    else:
        logger.info("start parallel processing with {} cores.".format(num_cores))
        pool = multiprocessing.Pool(processes=num_cores)

        demerr = np.zeros((net_obj.num_arcs, 1), dtype=np.float32)
        vel = np.zeros((net_obj.num_arcs, 1), dtype=np.float32)
        gamma = np.zeros((net_obj.num_arcs, 1), dtype=np.float32)

        num_cores = net_obj.num_arcs if num_cores > net_obj.num_arcs else num_cores  # avoids having more samples then
        # cores
        idx = ut.splitDatasetForParallelProcessing(num_samples=net_obj.num_arcs, num_cores=num_cores)

        args = [(
            idx_range,
            idx_range.shape[0],
            net_obj.phase[idx_range, :],
            net_obj.slant_range[idx_range],
            net_obj.loc_inc[idx_range],
            ifg_net_obj,
            wavelength,
            velocity_bound,
            demerr_bound,
            num_samples) for idx_range in idx]

        results = pool.map(func=launchAmbiguityFunctionSearch, iterable=args)

        # retrieve results
        for i, demerr_i, vel_i, gamma_i in results:
            demerr[i] = demerr_i
            vel[i] = vel_i
            gamma[i] = gamma_i

    m, s = divmod(time.time() - start_time, 60)
    logger.info("Finished temporal unwrapping.")
    logger.debug('time used: {:02.0f} mins {:02.1f} secs.'.format(m, s))
    return demerr, vel, gamma


def launchSpatialUnwrapping(parameters: tuple) -> tuple[np.ndarray, np.ndarray]:
    """LaunchSpatialUnwrapping.

    Parameters
    ----------
    parameters: tuple
        idx_range, num_ifgs, num_points, edges, phase

    Returns
    -------
    idx_range: np.ndarray
    unw_phase: np.ndarray
    """
    # Unpack the parameters
    (idx_range, num_ifgs, num_points, method, edges, phase) = parameters

    prog_bar = ptime.progressBar(maxValue=num_ifgs)

    unw_phase = np.zeros((num_points, num_ifgs), dtype=np.float32)

    # Perform the PUMA phase unwrapping
    for i in range(num_ifgs):
        if method == "puma":
            unw_phase[:, i] = unwrap_arbitrary(
                psi=phase[:, i],
                edges=edges,
                simplices=None,
                method="gc",
                period=2*np.pi,
                start_i=0,
                p=0.2
            )
        else:
            unw_phase[:, i] = unwrap_arbitrary(
                psi=phase[:, i],
                edges=edges,
                simplices=None,  # todo: compute simplices for ILP
                method="ilp",
                period=2*np.pi,
                start_i=0,
            )
        prog_bar.update(value=i + 1, every=1,
                        suffix='{}/{} ifgs unwrapped. '.format(i + 1, num_ifgs))

    unw_phase = unw_phase - np.mean(unw_phase, axis=0)
    return idx_range, unw_phase


def spatialUnwrapping(*, num_ifgs: int, num_points: int, phase: np.ndarray, edges: np.ndarray, method: str,
                      num_cores: int, logger: Logger):
    """Spatial unwrapping of interferograms for a set of points.

    Parameters
    ----------
    num_ifgs: int
        Number of interferograms.
    num_points: int
        Number of points.
    phase: np.ndarray
        Phase of the interferograms at the points.
    edges: np.ndarray
        Edges/arcs of the graph.
    method: str
        Method for spatial unwrapping (puma or ilp).
    num_cores: int
        Number of cores to be used in multiprocessing.
    logger: Logger
        Logging handler.

    Returns
    -------
    unw_phase: np.ndarray
        Unwrapped phase of the interferograms at the points.
    """
    msg = "#" * 10
    msg += f" SPATIAL UNWRAPPING: {method} "
    msg += "#" * 10
    logger.info(msg)

    start_time = time.time()

    if num_cores == 1:
        parameters = (
            np.arange(num_ifgs),
            num_ifgs,
            num_points,
            method,
            edges,
            phase
        )
        idx_range, unw_phase = launchSpatialUnwrapping(parameters=parameters)
    else:
        logger.info("start parallel processing with {} cores.".format(num_cores))
        pool = multiprocessing.Pool(processes=num_cores)

        unw_phase = np.zeros((num_points, num_ifgs), dtype=np.float32)
        num_cores = num_ifgs if num_cores > num_ifgs else num_cores
        # avoids having more samples than cores
        idx = ut.splitDatasetForParallelProcessing(num_samples=num_ifgs, num_cores=num_cores)

        args = [(
            idx_range,
            idx_range.shape[0],
            num_points,
            method,
            edges,
            phase[:, idx_range]) for idx_range in idx]
        results = pool.map(func=launchSpatialUnwrapping, iterable=args)

        # retrieve results
        for i, phase in results:
            unw_phase[:, i] = phase

    m, s = divmod(time.time() - start_time, 60)
    logger.debug('time used: {:02.0f} mins {:02.1f} secs.'.format(m, s))

    return unw_phase


def spatialParameterIntegrationIterative(*,
                                         val_arcs: np.ndarray,
                                         all_arcs: np.ndarray,
                                         coord_xy: np.ndarray,
                                         all_weights: np.ndarray,
                                         spatial_ref_idx: int = 0,
                                         res_tol: float = 1e-3,
                                         max_rm_fraction: float = 0.001,
                                         logger: Logger):
    """Unwrapping double-difference arc parameters spatially.

    The parameters at the arcs are integrated spatially to the points. The integration is done iteratively using
    least-squares by removing the arcs with the highest residuals in each iteration.
    The integration stops when the sum of the residuals is below a threshold.
    Function is adopted from StaMPS software (Hooper et al., 2007).

    Parameters
    ----------
    val_arcs: np.ndarray
        Value at the arcs (e.g. DEM error, velocity).
    all_arcs: np.ndarray
        Arcs of the spatial network.
    coord_xy: np.ndarray
        Radar coordinates of the points in the spatial network.
    all_weights: np.ndarray
        Weights of the arcs (e.g. temporal coherence from temporal unwrapping)
    spatial_ref_idx: int
        Index of the spatial reference point (default = 0). Can be arbitrary.
    res_tol: float
        Threshold on the sum of the residual phase (default = 1e-3). Convergence criterion.
    max_rm_fraction: float
        Fraction of the arcs that are removed in each iteration (default = 0.001).
    logger: Logger
        Logging handler

    Returns
    -------
    val_points: np.ndarray
        Estimated parameters at the points resulting from the integration of the parameters at the arcs.
    """
    all_arcs = np.array(all_arcs)
    num_points = coord_xy.shape[0]
    num_arcs = all_arcs.shape[0]

    # create design matrix
    a = np.zeros((num_arcs, num_points))
    for i in range(num_arcs):
        a[i, all_arcs[i][0]] = 1
        a[i, all_arcs[i][1]] = -1

    # find the number of arcs per point
    arcs_per_point = np.zeros(num_points, )

    for i in range(num_points):
        arcs_per_point[i] = np.where(a[:, i] != 0)[0].shape[0]

    # remove reference point from design matrix
    all_a = csr_matrix(all_weights * np.delete(a, spatial_ref_idx, 1))

    # don't even start if the network is not connected
    if structural_rank(all_a) < all_a.shape[1]:
        logger.exception("Spatial point network is not connected. Phase cannot be unwrapped!")
        raise Exception

    # set n_bad to maximum fraction of bad edges that can be removed
    n_bad = np.ceil(num_arcs * max_rm_fraction).astype(np.int16)

    # initialize output
    val_points = np.zeros((num_points,))
    points_idx = np.ones((num_points,), dtype=bool)
    points_idx[spatial_ref_idx] = False
    x_hat = np.zeros((num_points - 1,))

    start_time = time.time()

    arcs = all_arcs
    obv_vec = val_arcs.reshape(-1, ) * all_weights.reshape(-1, )
    a = all_a
    weights = all_weights
    num_arcs = obv_vec.size

    r = None
    num_arcs_save = None
    arcs_save = None
    a_save = None
    weights_save = None
    obv_vec_save = None
    i = 0
    while True:
        if structural_rank(a) >= a.shape[1]:
            x_hat[:] = lsqr(a, obv_vec)[0]

            # store the current version of variables, being able to go back to previous iteration if too many arcs
            # removed
            a_save = a
            obv_vec_save = obv_vec
            weights_save = weights
            arcs_save = arcs
            num_arcs_save = num_arcs

            # compute residuals
            r = obv_vec - np.matmul(a.toarray(), x_hat)

        else:  # network is not connected anymore, remove less psPoints and try again
            # x_hat = np.linalg.lstsq(a_save, obv_vec_save, rcond=None)[0]  # unclear: I think it is not necessary to
            # recompute the inversion.
            n_bad = np.ceil(n_bad / 10).astype(np.int16)  # remove less point

        if np.all(np.abs(r) < res_tol):
            break
        else:
            # drop arcs with the highest residuals, but only drop max one arc per point
            ps_w_dropped_arc = np.zeros((num_points,))
            good_arc_idx = np.ones((num_arcs_save,), dtype=bool)
            r_sort_idx = np.abs(r).argsort()[::-1]  # descending order, makes for loop easier

            for j in range(n_bad):  # remove arcs one by one
                bad_arc_idx = r_sort_idx[j]
                ps_idx0 = arcs_save[bad_arc_idx][0]
                ps_idx1 = arcs_save[bad_arc_idx][1]
                if (ps_w_dropped_arc[ps_idx0] == 0) and (ps_w_dropped_arc[
                                                             ps_idx1] == 0):  # if arc not already dropped for either
                    # point of current arc drop current arc
                    good_arc_idx[bad_arc_idx] = False
                    # mark both psPoints from the arc as having an arc dropped
                    ps_w_dropped_arc[ps_idx0] = 1
                    ps_w_dropped_arc[ps_idx1] = 1

            # update all variables for next iteration
            arcs = arcs_save[good_arc_idx, :]
            obv_vec = obv_vec_save[good_arc_idx]
            a = a_save[good_arc_idx, :]
            weights = weights_save[good_arc_idx]
            num_arcs = obv_vec.size

        i += 1

    val_points[points_idx] = x_hat

    m, s = divmod(time.time() - start_time, 60)
    logger.debug('time used: {:02.0f} mins {:02.1f} secs.\n'.format(m, s))

    return val_points


def spatialParameterIntegration(*,
                                val_arcs: np.ndarray,
                                arcs: np.ndarray,
                                coord_xy: np.ndarray,
                                weights: np.ndarray,
                                spatial_ref_idx: int = 0,
                                logger: Logger):
    """Unwrapping double-difference arc parameters spatially.

    The parameters at the arcs are integrated spatially to the points. The integration is done using least-squares.

    Parameters
    ----------
    val_arcs: np.ndarray
        Value at the arcs (e.g. DEM error, velocity).
    arcs: np.ndarray
        Arcs of the spatial network.
    coord_xy: np.ndarray
        Radar coordinates of the points in the spatial network.
    weights: np.ndarray
        Weights of the arcs (e.g. temporal coherence from temporal unwrapping)
    spatial_ref_idx: int
        Index of the spatial reference point (default = 0). Can be arbitrary.
    logger: Logger
        Logging handler

    Returns
    -------
    val_points: np.ndarray
        Estimated parameters at the points resulting from the integration of the parameters at the arcs.
    """
    arcs = np.array(arcs)
    num_points = coord_xy.shape[0]
    num_arcs = arcs.shape[0]

    # create design matrix
    design_mat = np.zeros((num_arcs, num_points))
    for i in range(num_arcs):
        design_mat[i, arcs[i][0]] = 1
        design_mat[i, arcs[i][1]] = -1

    # remove reference point from design matrix
    design_mat = csr_matrix(weights * np.delete(design_mat, spatial_ref_idx, 1))

    # don't even start if the network is not connected
    if structural_rank(design_mat) < design_mat.shape[1]:
        raise Exception("Spatial point network is not connected. Cannot integrate parameters spatially!")

    start_time = time.time()

    obv_vec = val_arcs.reshape(-1, ) * weights.reshape(-1, )

    x_hat = lsqr(design_mat, obv_vec)[0]

    m, s = divmod(time.time() - start_time, 60)
    logger.debug('time used: {:02.0f} mins {:02.1f} secs.'.format(m, s))

    val_points = np.zeros((num_points,))
    points_idx = np.ones((num_points,), dtype=bool)
    points_idx[spatial_ref_idx] = False
    val_points[points_idx] = x_hat

    return val_points


def computeNumArcsPerPoints(*, net_obj: Network, point_id: np.ndarray,
                            logger: Logger) -> tuple[np.ndarray, np.ndarray]:
    """Remove Points with less than specified number of arcs.

    Parameters
    ----------
    net_obj: Network
        The spatial Network object.
    point_id: np.ndarray
        ID of the points in the network.
    logger: Logger
        Logging handler.

    Returns
    -------
     design_mat: np.ndarray
        Design matrix of the spatial network
     arcs_per_point: np.ndarray
        Number of arcs that each point is connected with.
    """
    logger.info("Removal of arcs and PSC that cannot be tested.")

    num_points = point_id.shape[0]

    # create design matrix
    design_mat = np.zeros((net_obj.num_arcs, num_points))
    for i in range(net_obj.num_arcs):
        design_mat[i, net_obj.arcs[i][0]] = 1
        design_mat[i, net_obj.arcs[i][1]] = -1

    # find the number of arcs per point
    arcs_per_point = np.zeros(num_points, )

    for i in range(num_points):
        arcs_per_point[i] = np.where(design_mat[:, i] != 0)[0].shape[0]

    return design_mat, arcs_per_point


def computeAvgCoherencePerPoint(*, net_obj: Network, point_id: np.ndarray, logger: Logger) -> np.ndarray:
    """Compute the average coherence from all arcs that a point is connected with. Used to remove incoherent points.

    Parameters
    ----------
    net_obj: Network
        The Network object.
    point_id: np.ndarray
        ID of the points in the network.
    logger: Logger
        Logging handler.

    Returns
    -------
    mean_gamma_point: np.ndarray
        Average coherence per point
    """
    logger.info("Removal of points whose arcs are incoherent in average.")

    num_points = point_id.shape[0]

    # create design matrix
    a = np.zeros((net_obj.num_arcs, num_points))
    for i in range(net_obj.num_arcs):
        a[i, net_obj.arcs[i][0]] = net_obj.gamma[i]
        a[i, net_obj.arcs[i][1]] = net_obj.gamma[i]

    a[a == 0] = np.nan
    mean_gamma_point = np.nanmean(a, axis=0)

    return mean_gamma_point


def removeArcsByPointMask(*, net_obj: Union[Network, NetworkParameter], point_id: np.ndarray, coord_xy: np.ndarray,
                          p_mask: np.ndarray, design_mat: np.ndarray,
                          logger: Logger) -> tuple[Network, np.ndarray, np.ndarray, np.ndarray]:
    """Remove all entries related to the arc observations connected to the points which have a False value in p_mask.

    Parameters
    ----------
    net_obj: Network
        The Network object.
    point_id: np.ndarray
        ID of the points in the network.
    coord_xy: np.ndarray
        Radar coordinates of the points in the spatial network.
    p_mask: np.ndarray
        Boolean mask with True for points to keep, and False for points to remove.
    design_mat: np.ndarray
        Design matrix describing the relation between arcs and points.
    logger: Logger
        Logging handler.

    Returns
    -------
    net_obj: Network
        Network object without the removed arcs and points.
    point_id: np.ndarray
        ID of the points in the network without the removed points.
    coord_xy: np.ndarray
        Radar coordinates of the points in the spatial network without the removed points.
    design_mat: np.ndarray
        Design matrix describing the relation between arcs and points without the removed points and arcs.
    """
    # find respective arcs
    a_idx = list()
    for p_idx in np.where(~p_mask)[0]:
        a_idx.append(np.where(design_mat[:, p_idx] != 0)[0])

    if len(a_idx) != 0:
        a_idx = np.hstack(a_idx)
        a_mask = np.ones((net_obj.num_arcs,), dtype=np.bool_)
        a_mask[a_idx] = False
        net_obj.removeArcs(mask=a_mask)
        design_mat = design_mat[a_mask, :]
    else:
        a_idx = np.array(a_idx)  # so I can check the size

    # remove psPoints
    point_id = point_id[p_mask]
    design_mat = design_mat[:, p_mask]
    coord_xy = coord_xy[p_mask, :]

    # beside removing the arcs in "arcs", the tuple indices have to be changed to make them fit to new point indices
    for p_idx in np.sort(np.where(~p_mask)[0])[::-1]:
        net_obj.arcs[np.where((net_obj.arcs[:, 0] > p_idx)), 0] -= 1
        net_obj.arcs[np.where((net_obj.arcs[:, 1] > p_idx)), 1] -= 1

    logger.info("Removed {} arc(s) connected to the removed point(s)".format(a_idx.size))
    return net_obj, point_id, coord_xy, design_mat


def removeGrossOutliers(*, net_obj: Network, point_id: np.ndarray, coord_xy: np.ndarray, min_num_arc: int = 3,
                        quality_thrsh: float = 0.0,
                        logger: Logger) -> tuple[Network, np.ndarray, np.ndarray, np.ndarray]:
    """Remove both gross outliers which have many low quality arcs and points which are not well connected.

    Parameters
    ----------
    net_obj: Network
        The spatial Network object.
    point_id: np.ndarray
        ID of the points in the network.
    coord_xy: np.ndarray
        Radar coordinates of the points in the spatial network.
    min_num_arc: int
        Threshold on the minimal number of arcs per point. Default = 3.
    quality_thrsh: float
        Threshold on the temporal coherence of the arcs. Default = 0.0.
    logger: Logger
        Logging handler.

    Returns
    -------
    net_obj: Network
        Network object without the removed arcs and points.
    point_id: np.ndarray
        ID of the points in the network without the removed points.
    coord_xy: np.ndarray
        Radar coordinates of the points in the spatial network without the removed points.
    a: np.ndarray
        Design matrix describing the relation between arcs and points without the removed points and arcs.
    """
    logger.info("Detect points with low quality arcs (mean): < {}".format(quality_thrsh))
    mean_gamma_point = computeAvgCoherencePerPoint(net_obj=net_obj,
                                                   point_id=point_id, logger=logger)
    # not yet removed, because arcs are removed separately
    p_mask_mean_coh = (mean_gamma_point >= quality_thrsh).ravel()
    logger.info("Detected {} point(s) with mean coherence of all connected arcs < {} ".format(
        p_mask_mean_coh[~p_mask_mean_coh].shape[0], quality_thrsh))

    logger.info("Removal of low quality arcs: < {}".format(quality_thrsh))
    a_mask = (net_obj.gamma >= quality_thrsh).ravel()
    logger.info("Removed {} arc(s)".format(a_mask[~a_mask].shape[0]))
    net_obj.removeArcs(mask=a_mask)

    design_mat, arcs_per_point = computeNumArcsPerPoints(net_obj=net_obj, point_id=point_id, logger=logger)

    p_mask_num_arcs = (arcs_per_point >= min_num_arc).ravel()
    logger.info("Detected {} point(s) with less than {} arcs".format(p_mask_num_arcs[~p_mask_num_arcs].shape[0],
                                                                         min_num_arc))

    # remove them jointly
    p_mask = p_mask_num_arcs & p_mask_mean_coh
    logger.info("Remove {} point(s)".format(p_mask[~p_mask].shape[0]))
    net_obj, point_id, coord_xy, design_mat = removeArcsByPointMask(net_obj=net_obj, point_id=point_id,
                                                                    coord_xy=coord_xy, p_mask=p_mask,
                                                                    design_mat=design_mat, logger=logger)
    return net_obj, point_id, coord_xy, design_mat


def parameterBasedNoisyPointRemoval(*, net_par_obj: NetworkParameter, point_id: np.ndarray, coord_xy: np.ndarray,
                                    design_mat: np.ndarray, rmse_thrsh: float = 0.02, num_points_remove: int = 1,
                                    bmap_obj: AmplitudeImage = None, bool_plot: bool = False,
                                    logger: Logger):
    """Remove Points during spatial integration step if residuals at many connected arcs are high.

    The idea is similar to outlier removal in DePSI, but without hypothesis testing.
    It can be used as a preprocessing step to spatial integration.
    The points are removed based on the RMSE computed from the residuals of the parameters (DEM error, velocity) per
    arc. The point with the highest RMSE is removed in each iteration. The process stops when the maximum RMSE is below
    a threshold.


    Parameters
    ----------
    net_par_obj: NetworkParameter
        The spatial NetworkParameter object containing the parameters estimates at each arc.
    point_id: np.ndarray
        ID of the points in the network.
    coord_xy: np.ndarray
        Radar coordinates of the points in the spatial network.
    design_mat: np.ndarray
        Design matrix describing the relation between arcs and points.
    rmse_thrsh: float
        Threshold for the RMSE of the residuals per point. Default = 0.02.
    num_points_remove: int
        Number of points to remove in each iteration. Default = 1.
    bmap_obj: AmplitudeImage
        Basemap object for plotting. Default = None.
    bool_plot: bool
        Plot the RMSE per point. Default = False.
    logger: Logger
        Logging handler.

    Returns
    -------
    spatial_ref_id: int
        ID of the spatial reference point.
    point_id: np.ndarray
        ID of the points in the network without the removed points.
    net_par_obj: NetworkParameter
        The NetworkParameter object without the removed points.
    """
    msg = "#" * 10
    msg += " NOISY POINT REMOVAL BASED ON ARC PARAMETERS "
    msg += "#" * 10
    logger.info(msg)

    num_points = point_id.shape[0]

    logger.info("Selection of the reference PSC")
    # select one of the two pixels which are connected via the arc with the highest quality
    spatial_ref_idx = np.where(design_mat[np.argmax(net_par_obj.gamma), :] != 0)[0][0]
    coord_xy = np.delete(arr=coord_xy, obj=spatial_ref_idx, axis=0)
    spatial_ref_id = point_id[spatial_ref_idx]
    point_id = np.delete(arr=point_id, obj=spatial_ref_idx, axis=0)
    num_points -= 1

    # remove reference point from design matrix
    design_mat = net_par_obj.gamma * np.delete(arr=design_mat, obj=spatial_ref_idx, axis=1)

    logger.info("Spatial integration to detect noisy point")
    start_time = time.time()

    it_count = 0
    while True:
        logger.info("ITERATION: {}".format(it_count))
        design_mat = csr_matrix(design_mat)

        if structural_rank(design_mat) < design_mat.shape[1]:
            logger.error("Singular normal matrix. Network is no longer connected!")
            # point_id = np.sort(np.hstack([spatial_ref_id, point_id]))
            # return spatial_ref_id, point_id, net_par_obj
            raise ValueError
        # demerr
        obv_vec = net_par_obj.demerr.reshape(-1, )
        demerr_points = lsqr(design_mat.toarray(), obv_vec * net_par_obj.gamma.reshape(-1, ))[0]
        r_demerr = obv_vec - np.matmul(design_mat.toarray(), demerr_points)

        # vel
        obv_vec = net_par_obj.vel.reshape(-1, )
        vel_points = lsqr(design_mat.toarray(), obv_vec * net_par_obj.gamma.reshape(-1, ))[0]
        r_vel = obv_vec - np.matmul(design_mat.toarray(), vel_points)

        rmse_demerr = np.zeros((num_points,))
        rmse_vel = np.zeros((num_points,))
        for p in range(num_points):
            r_mask = design_mat[:, p].toarray() != 0
            rmse_demerr[p] = np.sqrt(np.mean(r_demerr[r_mask.ravel()].ravel() ** 2))
            rmse_vel[p] = np.sqrt(np.mean(r_vel[r_mask.ravel()].ravel() ** 2))

        rmse = rmse_vel.copy()
        max_rmse = np.max(rmse.ravel())
        logger.info("Maximum RMSE DEM correction: {:.2f} m".format(np.max(rmse_demerr.ravel())))
        logger.info("Maximum RMSE velocity: {:.4f} m / year".format(np.max(rmse_vel.ravel())))

        if bool_plot:
            # vel
            ax = bmap_obj.plot(logger=logger)
            sc = ax.scatter(coord_xy[:, 1], coord_xy[:, 0], c=rmse_vel * 1000, s=3.5,
                            cmap=plt.cm.get_cmap("autumn_r"), vmin=0, vmax=rmse_thrsh * 1000)
            plt.colorbar(sc, pad=0.03, shrink=0.5)
            ax.set_title("{}. iteration\nmean velocity - RMSE per point in [mm / year]".format(it_count))
            fig = ax.get_figure()
            plt.tight_layout()
            fig.savefig(join(dirname(net_par_obj.file_path), "pic", f"step_1_rmse_vel_{it_count}th_iter.png"),
                        dpi=300)
            plt.close(fig)

            # demerr
            ax = bmap_obj.plot(logger=logger)
            sc = ax.scatter(coord_xy[:, 1], coord_xy[:, 0], c=rmse_demerr, s=3.5,
                            cmap=plt.cm.get_cmap("autumn_r"))
            plt.colorbar(sc, pad=0.03, shrink=0.5)
            ax.set_title("{}. iteration\nDEM correction - RMSE per point in [m]".format(it_count))
            fig = ax.get_figure()
            plt.tight_layout()
            fig.savefig(join(dirname(net_par_obj.file_path), "pic",
                             f"step_1_rmse_dem_correction_{it_count}th_iter.png"),
                        dpi=300)
            plt.close(fig)

        if max_rmse <= rmse_thrsh:
            logger.info("No noisy pixels detected.")
            break

        # remove point with highest rmse
        p_mask = np.ones((num_points,), dtype=np.bool_)
        p_mask[np.argsort(rmse)[::-1][:num_points_remove]] = False  # see description of function removeArcsByPointMask
        net_par_obj, point_id, coord_xy, design_mat = removeArcsByPointMask(net_obj=net_par_obj, point_id=point_id,
                                                                            coord_xy=coord_xy, p_mask=p_mask,
                                                                            design_mat=design_mat.toarray(),
                                                                            logger=logger)
        num_points -= num_points_remove
        it_count += 1

    m, s = divmod(time.time() - start_time, 60)
    logger.debug('time used: {:02.0f} mins {:02.1f} secs.'.format(m, s))

    # add spatialRefIdx back to point_id
    point_id = np.sort(np.hstack([spatial_ref_id, point_id]))
    return spatial_ref_id, point_id, net_par_obj
