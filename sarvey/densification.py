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

"""Densification module for SARvey."""
import time
import multiprocessing
import numpy as np
from scipy.spatial import KDTree
from logging import Logger

from mintpy.utils import ptime

from sarvey.unwrapping import oneDimSearchTemporalCoherence
from sarvey.objects import Points
import sarvey.utils as ut


def densificationInitializer(tree_p1: KDTree, point2_obj: Points, demod_phase1: np.ndarray):
    """DensificationInitializer.

    Sets values to global variables for parallel processing.

    Parameters
    ----------
    tree_p1 : KDTree
        KDTree of the first-order network
    point2_obj : Points
        Points object with second-order points
    demod_phase1 : np.ndarray
        demodulated phase of the first-order network
    """
    global global_tree_p1
    global global_point2_obj
    global global_demod_phase1

    global_tree_p1 = tree_p1
    global_point2_obj = point2_obj
    global_demod_phase1 = demod_phase1


def launchDensifyNetworkConsistencyCheck(args: tuple):
    """LaunchDensifyNetworkConsistencyCheck.

    Launches the densification of the network with second-order points inside parallel processing.

    Parameters
    ----------
    args : tuple
        Tuple with the following parameters:

        idx_range : np.ndarray
            Array with the indices of the second-order points
        num_points : int
            Number of second-order points
        num_conn_p1 : int
            Number of nearest points in the first-order network
        max_dist_p1 : float
            Maximum allowed distance to the nearest points in the first-order network
        velocity_bound : float
            Bound for the velocity estimate in temporal unwrapping
        demerr_bound : float
            Bound for the DEM error estimate in temporal unwrapping
        num_samples : int
            Number of samples for the search of the optimal parameters

    Returns
    -------
    idx_range : np.ndarray
        Array with the indices of the second-order points
    demerr_p2 : np.ndarray
        DEM error array of the second-order points
    vel_p2 : np.ndarray
        Velocity array of the second-order points
    gamma_p2 : np.ndarray
        Estimated temporal coherence array of the second-order points resulting from temporal unwrapping
    """
    (idx_range, num_points, num_conn_p1, max_dist_p1, velocity_bound, demerr_bound, num_samples) = args

    counter = 0
    prog_bar = ptime.progressBar(maxValue=num_points)

    # initialize output
    demerr_p2 = np.zeros((num_points,), dtype=np.float32)
    vel_p2 = np.zeros((num_points,), dtype=np.float32)
    gamma_p2 = np.zeros((num_points,), dtype=np.float32)

    design_mat = np.zeros((global_point2_obj.ifg_net_obj.num_ifgs, 2), dtype=np.float32)

    demerr_range = np.linspace(-demerr_bound, demerr_bound, num_samples)
    vel_range = np.linspace(-velocity_bound, velocity_bound, num_samples)

    factor = 4 * np.pi / global_point2_obj.wavelength

    for idx in range(num_points):
        p2 = idx_range[idx]
        # nearest points in p1
        dist, nearest_p1 = global_tree_p1.query([global_point2_obj.coord_utm[p2, 0],
                                                 global_point2_obj.coord_utm[p2, 1]], k=num_conn_p1)
        mask = (dist < max_dist_p1) & (dist != 0)
        mask[:3] = True  # ensure that always at least the three closest points are used
        nearest_p1 = nearest_p1[mask]

        # compute arc observations to nearest points
        arc_phase_p1 = np.angle(np.exp(1j * global_point2_obj.phase[p2, :]) *
                                np.conjugate(np.exp(1j * global_demod_phase1[nearest_p1, :])))

        design_mat[:, 0] = (factor * global_point2_obj.ifg_net_obj.pbase_ifg
                            / (global_point2_obj.slant_range[p2] * np.sin(global_point2_obj.loc_inc[p2])))
        design_mat[:, 1] = factor * global_point2_obj.ifg_net_obj.tbase_ifg

        demerr_p2[idx], vel_p2[idx], gamma_p2[idx] = oneDimSearchTemporalCoherence(
            demerr_range=demerr_range,
            vel_range=vel_range,
            obs_phase=arc_phase_p1,
            design_mat=design_mat
        )

        prog_bar.update(counter + 1, every=np.int16(200),
                        suffix='{}/{} points'.format(counter + 1, num_points))
        counter += 1

    return idx_range, demerr_p2, vel_p2, gamma_p2


def densifyNetwork(*, point1_obj: Points, vel_p1: np.ndarray, demerr_p1: np.ndarray, point2_obj: Points,
                   num_conn_p1: int, max_dist_p1: float, velocity_bound: float, demerr_bound: float,
                   num_samples: int, num_cores: int = 1, logger: Logger):
    """DensifyNetwork.

    Densifies the network with second-order points by connecting the second-order points to the closest points in the
    first-order network.

    Parameters
    ----------
    point1_obj : Points
        Points object with first-order points
    vel_p1 : np.ndarray
        Velocity array of the first-order points
    demerr_p1 : np.ndarray
        DEM error array of the first-order points
    point2_obj : Points
        Points object with second-order points
    num_conn_p1 : int
        Number of nearest points in the first-order network
    max_dist_p1 : float
        Maximum allowed distance to the nearest points in the first-order network
    velocity_bound : float
        Bound for the velocity estimate in temporal unwrapping
    demerr_bound : float
        Bound for the DEM error estimate in temporal unwrapping
    num_samples : int
        Number of samples for the search of the optimal parameters
    num_cores : int
        Number of cores for parallel processing (default: 1)
    logger : Logger
        Logger object

    Returns
    -------
    demerr_p2 : np.ndarray
        DEM error array of the second-order points
    vel_p2 : np.ndarray
        Velocity array of the second-order points
    gamma_p2 : np.ndarray
        Estimated temporal coherence array of the second-order points resulting from temporal unwrapping
    """
    msg = "#" * 10
    msg += " DENSIFICATION WITH SECOND-ORDER POINTS "
    msg += "#" * 10
    logger.info(msg=msg)
    start_time = time.time()

    # find the closest points from first-order network
    tree_p1 = KDTree(data=point1_obj.coord_utm)

    # remove parameters from wrapped phase
    pred_phase_demerr, pred_phase_vel = ut.predictPhase(
        obj=point1_obj,
        vel=vel_p1, demerr=demerr_p1,
        ifg_space=True, logger=logger
    )
    pred_phase = pred_phase_demerr + pred_phase_vel

    # Note: for small baselines it does not make a difference if re-wrapping the phase difference or not.
    # However, for long baselines (like in the star network) it does make a difference. Leijen (2014) does not re-wrap
    # the arc double differences to be able to test the ambiguities. Kampes (2006) does re-wrap, but is testing based
    # on the estimated parameters. Hence, it doesn't make a difference for him. Not re-wrapping can be a starting point
    # for triangle-based temporal unwrapping.
    # demod_phase1 = np.angle(np.exp(1j * point1_obj.phase) * np.conjugate(np.exp(1j * pred_phase)))  # re-wrapping
    demod_phase1 = point1_obj.phase - pred_phase  # not re-wrapping

    # initialize output
    init_args = (tree_p1, point2_obj, demod_phase1)

    if num_cores == 1:
        densificationInitializer(tree_p1=tree_p1, point2_obj=point2_obj, demod_phase1=demod_phase1)
        args = (np.arange(point2_obj.num_points), point2_obj.num_points, num_conn_p1, max_dist_p1,
                velocity_bound, demerr_bound, num_samples)
        idx_range, demerr_p2, vel_p2, gamma_p2 = launchDensifyNetworkConsistencyCheck(args)
    else:
        with multiprocessing.Pool(num_cores, initializer=densificationInitializer, initargs=init_args) as pool:
            logger.info(msg="start parallel processing with {} cores.".format(num_cores))
            num_cores = point2_obj.num_points if num_cores > point2_obj.num_points else num_cores
            # avoids having less samples than cores
            idx = ut.splitDatasetForParallelProcessing(num_samples=point2_obj.num_points, num_cores=num_cores)
            args = [(
                idx_range,
                idx_range.shape[0],
                num_conn_p1,
                max_dist_p1,
                velocity_bound,
                demerr_bound,
                num_samples
            ) for idx_range in idx]

            results = pool.map_async(launchDensifyNetworkConsistencyCheck, args, chunksize=1)
            while True:
                time.sleep(5)
                if results.ready():
                    results = results.get()
                    break
            # needed to make coverage work in multiprocessing (not sure what that means. copied from package Arosics).
            pool.close()
            pool.join()

        demerr_p2 = np.zeros((point2_obj.num_points,), dtype=np.float32)
        vel_p2 = np.zeros((point2_obj.num_points,), dtype=np.float32)
        gamma_p2 = np.zeros((point2_obj.num_points,), dtype=np.float32)

        # retrieve results
        for i, demerr_i, vel_i, gamma_i in results:
            demerr_p2[i] = demerr_i
            vel_p2[i] = vel_i
            gamma_p2[i] = gamma_i

    m, s = divmod(time.time() - start_time, 60)
    logger.debug(msg='time used: {:02.0f} mins {:02.1f} secs.\n'.format(m, s))

    # combine p1 and p2 parameters and bring them in correct order using point_id
    sort_idx = np.argsort(np.append(point1_obj.point_id, point2_obj.point_id))
    demerr_p2 = np.append(demerr_p1, demerr_p2)  # add gamma=1 for p1 pixels
    vel_p2 = np.append(vel_p1, vel_p2)
    gamma_p2 = np.append(np.ones_like(point1_obj.point_id), gamma_p2)  # add gamma=1 for p1 pixels

    demerr_p2 = demerr_p2[sort_idx]
    vel_p2 = vel_p2[sort_idx]
    gamma_p2 = gamma_p2[sort_idx]
    return demerr_p2, vel_p2, gamma_p2
