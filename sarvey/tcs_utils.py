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
"""Utils module for TCS processing with SARvey."""
import logging
import time
import numpy as np
from logging import Logger
import multiprocessing

from mintpy.utils import ptime

from sarvey.objects import Points
from sarvey.triangulation import PointNetworkTriangulation
from sarvey import utils as ut


def launchPreTriangulateEachIfg(parameters: tuple) -> tuple[np.ndarray, list[np.ndarray]]:
    """LaunchPreTriangulateEachIfg.

    Parameters
    ----------
    parameters: tuple
        idx_range, num_ifgs, coord_xy, lifetime_ifgs

    Returns
    -------
    idx_range: np.ndarray
        Index range.
    edges_per_ifg: list[np.ndarray]
        Edges per interferogram resulting from triangulation.
    """
    # Unpack the parameters
    (idx_range, num_ifgs, coord_xy, lifetime_ifgs) = parameters

    edges_per_ifg = list()
    prog_bar = ptime.progressBar(maxValue=num_ifgs)

    for idx in range(num_ifgs):
        triang_obj = PointNetworkTriangulation(
            coord_xy=coord_xy[lifetime_ifgs[:, idx]],
            coord_map_xy=None,
            logger=logging.Logger("unused"),
            verbose=False
        )
        triang_obj.triangulateGlobal()  # if coord_map is not given, only global delaunay and knn can be calculated
        edges_per_ifg.append(triang_obj.getArcsFromAdjMat())

        prog_bar.update(value=idx + 1, every=5, suffix=f'{idx + 1}/{num_ifgs} ifgs triangulated')

    return idx_range, edges_per_ifg


def preTriangulateEachIfg(*, coord_xy: np.ndarray, lifetime_ifgs: np.ndarray, num_cores: int, logger: Logger):
    """PreTriangulateEachIfg.

    Apply Delaunay triangulation to create arcs for the spatial phase unwrapping. This is done for each interferogram
    separately as different sets of points are coherent in each interferogram.

    Parameters
    ----------
    coord_xy: np.ndarray
        XY coordinates of the points that shall be triangulated.
    lifetime_ifgs: np.ndarray
        Indicates coherent lifetime of points per interferograms (num_points x num_ifgs).
    num_cores: int
        Number of cores to use for parallel processing.
    logger: Logger
        Logger object.

    Returns
    -------
    edges_per_ifg: list[np.ndarray]
        Edges per interferogram resulting from triangulation (list[(num_arcs x 2)]).
    """
    edges_per_ifg = list()
    num_ifgs = lifetime_ifgs.shape[1]

    start_time = time.time()

    if num_cores == 1:
        parameters = (
            np.arange(num_ifgs),
            num_ifgs,
            coord_xy,
            lifetime_ifgs
        )
        idx_range, edges_per_ifg = launchPreTriangulateEachIfg(parameters=parameters)
    else:
        logger.info(msg="start parallel processing with {} cores.".format(num_cores))
        pool = multiprocessing.Pool(processes=num_cores)

        num_cores = num_ifgs if num_cores > num_ifgs else num_cores  # avoids having more samples than cores
        idx = ut.splitDatasetForParallelProcessing(num_samples=num_ifgs, num_cores=num_cores)

        args = [(
            idx_range,
            idx_range.shape[0],
            coord_xy,
            lifetime_ifgs[:, idx_range]) for idx_range in idx]
        results = pool.map(func=launchPreTriangulateEachIfg, iterable=args)

        # retrieve results
        results = sorted(results, key=lambda x: x[0][0])  # ensure the correct order of the results
        for i, edges in results:
            edges_per_ifg.extend(edges)

    m, s = divmod(time.time() - start_time, 60)
    logger.debug(msg='time used: {:02.0f} mins {:02.1f} secs.'.format(m, s))
    return edges_per_ifg


def selectCoherentImagesAndIfgs(*, change_index_map: np.ndarray, subset_index_map: np.ndarray, point_obj: Points,
                                mask_ccs: np.ndarray, mask_tcs: np.ndarray):
    """Get masks for coherent lifetime of TCS."""
    mask_p2 = point_obj.createMask()

    ccs_tcs_map = np.zeros_like(mask_p2, dtype=np.int8)  # 0: not selected
    ccs_tcs_map[mask_tcs] = 1  # 1: TCS
    ccs_tcs_map[mask_ccs] = 2  # 2: CCS

    ccs_tcs_points = ccs_tcs_map[mask_p2]

    lifetime_images = np.zeros((point_obj.num_points, point_obj.ifg_net_obj.num_images), dtype=np.bool_)
    lifetime_images[ccs_tcs_points == 2, :] = True  # CCS are coherent the whole time

    # add individual coherent lifetime of each TCS
    subset_index = subset_index_map[mask_tcs]
    change_idx = change_index_map[mask_tcs]
    tcs_point_idx = np.where(ccs_tcs_points == 1)[0]  # index of TCS points within the data of point_obj

    for i in range(mask_tcs.sum()):
        if subset_index[i] == 0:  # first part of time span is coherent
            lifetime_images[tcs_point_idx[i], :change_idx[i] + 1] = True
        else:  # second part of time span is coherent
            lifetime_images[tcs_point_idx[i], change_idx[i] + 1:] = True

    # convert coh_lifetime from images to ifgs
    design_mat = point_obj.ifg_net_obj.getDesignMatrix()
    """idea: for an interferogram to be within the coherent lifetime, it has to have two entries inside the row (ifg) of
    the design matrix for columns that belong to the coherent part (lifetime_images)"""
    lifetime_ifgs = np.zeros((point_obj.num_points, point_obj.ifg_net_obj.num_ifgs), dtype=np.bool_)
    lifetime_ifgs[ccs_tcs_points == 2, :] = True  # CCS are coherent the whole time
    for i in range(mask_tcs.sum()):
        num_imgs_per_ifg = np.sum(np.abs(design_mat[:, lifetime_images[tcs_point_idx[i], :]]), axis=1)
        lifetime_ifgs[tcs_point_idx[i], num_imgs_per_ifg == 2] = True

    return lifetime_images, lifetime_ifgs
