#!/usr/bin/env python

# SARvey - A multitemporal InSAR time series tool for the derivation of displacements.
#
# Copyright (C) 2021-2026 Andreas Piter (IPI Hannover, piter@ipi.uni-hannover.de)
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
import time
import networkx as nx
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
from sarvey.objects import Network, NetworkParameter


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


def temporalUnwrapping(*, ifg_net_obj: IfgNetwork, net_obj: Network, wavelength: float, velocity_bound: float,
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
    logger.info(msg=msg)

    start_time = time.time()

    if num_cores == 1:
        args = (
            np.arange(net_obj.num_arcs), net_obj.num_arcs, net_obj.phase,
            net_obj.slant_range, net_obj.loc_inc, ifg_net_obj, wavelength, velocity_bound, demerr_bound, num_samples)
        arc_idx_range, demerr, vel, gamma = launchAmbiguityFunctionSearch(parameters=args)
    else:
        logger.info(msg="start parallel processing with {} cores.".format(num_cores))

        demerr = np.zeros((net_obj.num_arcs, 1), dtype=np.float32)
        vel = np.zeros((net_obj.num_arcs, 1), dtype=np.float32)
        gamma = np.zeros((net_obj.num_arcs, 1), dtype=np.float32)

        num_cores = net_obj.num_arcs if num_cores > net_obj.num_arcs else num_cores  # avoids having more samples
        # then cores
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

        with multiprocessing.Pool(processes=num_cores) as pool:
            results = pool.map(func=launchAmbiguityFunctionSearch, iterable=args)

        # retrieve results
        for i, demerr_i, vel_i, gamma_i in results:
            demerr[i] = demerr_i
            vel[i] = vel_i
            gamma[i] = gamma_i

    m, s = divmod(time.time() - start_time, 60)
    logger.info(msg="Finished temporal unwrapping.")
    logger.debug(msg='time used: {:02.0f} mins {:02.1f} secs.'.format(m, s))
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
    logger.info(msg=msg)

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
        logger.info(msg="start parallel processing with {} cores.".format(num_cores))

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

        with multiprocessing.Pool(processes=num_cores) as pool:
            results = pool.map(func=launchSpatialUnwrapping, iterable=args)

        # retrieve results
        for i, phase in results:
            unw_phase[:, i] = phase

    m, s = divmod(time.time() - start_time, 60)
    logger.debug(msg='time used: {:02.0f} mins {:02.1f} secs.'.format(m, s))

    return unw_phase


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
    logger.debug(msg='time used: {:02.0f} mins {:02.1f} secs.'.format(m, s))

    val_points = np.zeros((num_points,))
    points_idx = np.ones((num_points,), dtype=bool)
    points_idx[spatial_ref_idx] = False
    val_points[points_idx] = x_hat

    return val_points


def removeBadPointsIteratively(*, net_obj: NetworkParameter, point_id: np.ndarray,
                               quality_thrsh: float, logger: Logger) -> [NetworkParameter, np.ndarray]:
    """
    Remove bad points from a network. Points with many low-quality arcs are removed iteratively.

    Parameters
    ----------
    net_obj: NetworkParameter
        The NetworkParameter object.
    point_id: np.ndarray
        ID of the points in the network.
    quality_thrsh: float
        Threshold on the temporal coherence of the arcs (edge weights).
    logger: Logger
        Logging handler.

    Returns
    -------
    net_obj: NetworkParameter
        NetworkParameter object without the removed points and arcs.
    point_id: np.ndarray
        ID of the points in the network after the removal of bad points.
    """
    logger.info(msg="Remove points with arcs that have a median temporal coherence < {}".format(quality_thrsh))

    graph = nx.DiGraph()
    graph.add_nodes_from(
        [(i, {'point_id': id}) for (i, id) in enumerate(point_id)]
    )
    graph.add_edges_from(
        [(arc[0], arc[1], {'weight': net_obj.gamma[idx], 'arc_idx': idx}) for idx, arc in enumerate(net_obj.arcs)]
    )

    median_coherence = {
        u: np.nanmedian([graph[u][v]['weight'] for v in graph.successors(u)] +
                        [graph[v][u]['weight'] for v in graph.predecessors(u)])
        for u in graph.nodes()
    }

    while True:
        worst_node = min(median_coherence, key=median_coherence.get)

        if median_coherence[worst_node] >= quality_thrsh:
            break

        affected_nodes = set(graph.successors(worst_node)) | set(graph.predecessors(worst_node))
        for u in affected_nodes:
            median_coherence[u] = np.nanmedian([graph[u][v]['weight'] for v in graph.successors(u)] +
                                               [graph[v][u]['weight'] for v in graph.predecessors(u)])

        graph.remove_node(worst_node)
        logger.debug("Removing point %d with median coherence %.2f",
                     point_id[worst_node], median_coherence[worst_node])

        del median_coherence[worst_node]

    lookup_dict = {node: index for index, node in enumerate(graph.nodes)}
    new_arc_list = [(lookup_dict[edge[0]], lookup_dict[edge[1]]) for edge in graph.edges]
    new_point_id = [graph.nodes[node]['point_id'] for node in graph.nodes()]

    logger.debug("Number of nodes after/before removal due to low temporal coherence: %d / %d",
                 len(new_point_id), len(point_id))
    logger.debug("Number of points removed due to low temporal coherence: %d", len(point_id) - len(new_point_id))

    arc_idx = [graph.edges[edge]['arc_idx'] for edge in graph.edges()]
    net_obj.arcs = np.array(new_arc_list, dtype=np.int64)
    net_obj.gamma = net_obj.gamma[arc_idx]
    net_obj.vel = net_obj.vel[arc_idx]
    net_obj.demerr = net_obj.demerr[arc_idx]
    net_obj.loc_inc = net_obj.loc_inc[arc_idx]
    net_obj.slant_range = net_obj.slant_range[arc_idx]
    net_obj.phase = net_obj.phase[arc_idx, :]
    net_obj.num_arcs = len(new_arc_list)
    point_id = new_point_id

    # log values after bad point removal
    logger.debug("[Min, Max] temporal coherence of points after bad point removal: [%.3f, %.3f]",
                 np.min(net_obj.gamma), np.max(net_obj.gamma))
    logger.debug("[Min, Max] velocity of points after bad point removal: [%.3f, %.3f]",
                 np.min(net_obj.vel), np.max(net_obj.vel))
    logger.debug("[Min, Max] DEM residual of points after bad point removal: [%.3f, %.3f]",
                 np.min(net_obj.demerr), np.max(net_obj.demerr))
    logger.debug("[Min, Max] incidence angle of points after bad point removal: [%.3f, %.3f]",
                 np.min(net_obj.loc_inc), np.max(net_obj.loc_inc))
    logger.debug("[Min, Max] slant range of points after bad point removal: [%.3f, %.3f]",
                 np.min(net_obj.slant_range), np.max(net_obj.slant_range))
    logger.debug("[Min, Max] phase of points after bad point removal: [%.3f, %.3f]",
                 np.min(net_obj.phase), np.max(net_obj.phase))

    logger.info(msg="Finished removing bad points.")
    return net_obj, point_id


def removeBadArcsIteratively(*,
                             net_obj: NetworkParameter,
                             quality_thrsh: float = 0.0,
                             logger: Logger) -> NetworkParameter:
    """Remove bad arcs iteratively from network based on quality threshold, preserving the minimum spanning tree.

    Parameters
    ----------
    net_obj: NetworkParameter
        The spatial NetworkParameter object.
    quality_thrsh: float
        Threshold on the temporal coherence of the arcs. Default = 0.0.
    logger: Logger
        Logging handler.

    Returns
    -------
    net_obj: NetworkParameter
        NetworkParameter object without the removed arcs.
    """
    logger.info(msg="Iteratively removing bad arcs with quality < {}".format(quality_thrsh))

    graph = nx.Graph()
    for idx, arc in enumerate(net_obj.arcs):
        graph.add_edge(arc[0], arc[1], weight=1-net_obj.gamma[idx])

    mst = nx.minimum_spanning_tree(graph, algorithm="kruskal")
    mst_edges = set(mst.edges)

    # Identify bad arcs that are not part of the MST
    bad_arc_mask = (net_obj.gamma < quality_thrsh).ravel()
    bad_arcs = [
        (arc[0], arc[1]) for idx, arc in enumerate(net_obj.arcs)
        if bad_arc_mask[idx] and (arc[0], arc[1]) not in mst_edges and (arc[1], arc[0]) not in mst_edges
    ]
    logger.info(msg="Removing {} bad arc(s)".format(len(bad_arcs)))

    # Remove the bad arcs
    bad_arc_indices = [
        idx for idx, arc in enumerate(net_obj.arcs)
        if (arc[0], arc[1]) in bad_arcs or (arc[1], arc[0]) in bad_arcs
    ]
    mask = np.ones(net_obj.num_arcs, dtype=bool)
    mask[bad_arc_indices] = False
    net_obj.removeArcs(mask=mask)

    return net_obj
