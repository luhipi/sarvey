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
from os.path import join
import time
import multiprocessing
import numpy as np
from scipy.spatial import KDTree
from logging import Logger
import matplotlib.pyplot as plt

from mintpy.utils import ptime, readfile
from miaplpy.objects.slcStack import slcStack

from sarvey.unwrapping import oneDimSearchTemporalCoherence
from sarvey.objects import Points, AmplitudeImage, BaseStack, ApsParameters
import sarvey.utils as ut
from sarvey.preparation import selectPixels, createTimeMaskFromDates
from sarvey.filtering import applySimpleInterpolationToP2, applySpatialFilteringToP2


def densificationInitializer(tree_p1: KDTree, tree_p2: KDTree, point2_obj: Points, demod_phase1: np.ndarray):
    """DensificationInitializer.

    Sets values to global variables for parallel processing.

    Parameters
    ----------
    tree_p1 : KDTree
        KDTree of the first-order network
    tree_p2 : KDTree
        KDTree of the second-order network
    point2_obj : Points
        Points object with second-order points
    demod_phase1 : np.ndarray
        demodulated phase of the first-order network
    """
    global global_tree_p1
    global global_tree_p2
    global global_point2_obj
    global global_demod_phase1

    global_tree_p1 = tree_p1
    global_tree_p2 = tree_p2
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
        num_conn_p2 : int
            Number of nearest points in the second-order network
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
    mean_gamma : np.ndarray
        Mean temporal coherence array of the second-order points resulting from consistency check with neighbourhood
    """
    (idx_range, num_points, num_conn_p1, num_conn_p2, max_dist_p1, velocity_bound, demerr_bound, num_samples) = args

    counter = 0
    prog_bar = ptime.progressBar(maxValue=num_points)

    # initialize output
    demerr_p2 = np.zeros((num_points,), dtype=np.float32)
    vel_p2 = np.zeros((num_points,), dtype=np.float32)
    gamma_p2 = np.zeros((num_points,), dtype=np.float32)
    mean_gamma = np.zeros((num_points,), dtype=np.float32)

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

        if num_conn_p2 > 0:
            # nearest points in p2 (for consistency check with neighbourhood)
            dist, nearest_p2 = global_tree_p2.query([global_point2_obj.coord_utm[p2, 0],
                                                     global_point2_obj.coord_utm[p2, 1]], k=num_conn_p2)
            mask = (dist != 0)  # remove the point itself
            nearest_p2 = nearest_p2[mask]
            num_n2 = mask[mask].shape[0]

            pred_phase_p2 = ut.predictPhaseSingle(
                vel=vel_p2[idx], demerr=demerr_p2[idx],
                slant_range=global_point2_obj.slant_range[p2], loc_inc=global_point2_obj.loc_inc[p2],
                ifg_net_obj=global_point2_obj.ifg_net_obj, wavelength=global_point2_obj.wavelength,
                only_vel=False, ifg_space=True
            )
            demod_phase2 = np.angle(np.exp(1j * global_point2_obj.phase[p2]) * np.conjugate(np.exp(1j * pred_phase_p2)))

            gamma_n2 = np.zeros((num_n2,), dtype=np.float32)
            for n_idx, n2 in enumerate(nearest_p2):
                arc_phase_p2 = np.angle(np.exp(1j * demod_phase2) *
                                        np.conjugate(np.exp(1j * global_point2_obj.phase[n2, :])))

                gamma_n2[n_idx] = oneDimSearchTemporalCoherence(
                    demerr_range=demerr_range,
                    vel_range=vel_range,
                    obs_phase=arc_phase_p2,
                    design_mat=design_mat
                )[-1]

            mean_gamma[idx] = np.average(gamma_n2, weights=1/(dist[mask] ** 2))

        prog_bar.update(counter + 1, every=np.int16(200),
                        suffix='{}/{} points'.format(counter + 1, num_points))
        counter += 1

    return idx_range, demerr_p2, vel_p2, gamma_p2, mean_gamma


def densifyNetwork(*, point1_obj: Points, vel_p1: np.ndarray, demerr_p1: np.ndarray, point2_obj: Points,
                   num_conn_p1: int, num_conn_p2: int, max_dist_p1: float, velocity_bound: float, demerr_bound: float,
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
    num_conn_p2 : int
        Number of nearest points in the second-order network
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
    mean_gamma : np.ndarray
        Mean temporal coherence array of the second-order points resulting from consistency check with neighbourhood

    """
    msg = "#" * 10
    msg += " DENSIFICATION WITH SECOND-ORDER POINTS "
    msg += "#" * 10
    logger.info(msg=msg)
    start_time = time.time()

    # find the closest points from first-order network
    tree_p1 = KDTree(data=point1_obj.coord_utm)
    tree_p2 = KDTree(data=point2_obj.coord_utm)

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
    init_args = (tree_p1, tree_p2, point2_obj, demod_phase1)

    if num_cores == 1:
        densificationInitializer(tree_p1=tree_p1, tree_p2=tree_p2, point2_obj=point2_obj, demod_phase1=demod_phase1)
        args = (np.arange(point2_obj.num_points), point2_obj.num_points, num_conn_p1, num_conn_p2, max_dist_p1,
                velocity_bound, demerr_bound, num_samples)
        idx_range, demerr_p2, vel_p2, gamma_p2, mean_gamma = launchDensifyNetworkConsistencyCheck(args)
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
                num_conn_p2,
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
        mean_gamma = np.zeros((point2_obj.num_points,), dtype=np.float32)

        # retrieve results
        for i, demerr_i, vel_i, gamma_i, mean_gamma_i in results:
            demerr_p2[i] = demerr_i
            vel_p2[i] = vel_i
            gamma_p2[i] = gamma_i
            mean_gamma[i] = mean_gamma_i

    m, s = divmod(time.time() - start_time, 60)
    logger.debug(msg='time used: {:02.0f} mins {:02.1f} secs.\n'.format(m, s))

    # combine p1 and p2 parameters and bring them in correct order using point_id
    sort_idx = np.argsort(np.append(point1_obj.point_id, point2_obj.point_id))
    demerr_p2 = np.append(demerr_p1, demerr_p2)  # add gamma=1 for p1 pixels
    vel_p2 = np.append(vel_p1, vel_p2)
    gamma_p2 = np.append(np.ones_like(point1_obj.point_id), gamma_p2)  # add gamma=1 for p1 pixels
    mean_gamma = np.append(np.ones_like(point1_obj.point_id), mean_gamma)  # add mean_gamma=1 for p1 pixels

    demerr_p2 = demerr_p2[sort_idx]
    vel_p2 = vel_p2[sort_idx]
    gamma_p2 = gamma_p2[sort_idx]
    mean_gamma = mean_gamma[sort_idx]
    return demerr_p2, vel_p2, gamma_p2, mean_gamma


def selectP2(*, output_path: str, config: dict, logger: Logger):
    """Select second order points and interpolate APS to them.

    Parameters
    ----------
    output_path : str
        path to the directory with processing files.
    config : dict
        dictionay containing parameters.
    logger : Logger
        Logger instance for logging messages.

    Returns
    -------
    point2_obj : Points
        Points object with second-order points
    aps2_obj: Points
        Points object with APS estimation for second-order points
    """
    # this function is mainly adapted from sarvey.processing.runFiltering
    # TODO: Directly pass input parameters instead of config dictionary
    coherence_p2 = config.densification.coherence_p2
    logger.debug("Start Selecting 2nd order points.")

    coh_value = int(coherence_p2 * 100)

    p1_file = join(output_path, "p1_ts_filt.h5")
    logger.debug(f"Create a mask based on points in file {p1_file}.")
    point1_obj = Points(file_path=p1_file, logger=logger)
    point1_obj.open(path_inputs=config.data_directories.path_inputs)
    p1_mask = point1_obj.createMask()  # used later for selecting psCand2 when a spatial mask AOI is given.
    logger.debug(f"Mask created with {p1_mask.sum()} selected points.")

    bmap_file = join(output_path, "background_map.h5")
    logger.debug(f"Use file {bmap_file} to create a mask invalid points.")
    bmap_obj = AmplitudeImage(file_path=bmap_file)

    point_id_img = np.arange(0, point1_obj.length * point1_obj.width).reshape(
        (point1_obj.length, point1_obj.width))

    logger.debug(f"Select second-order points based on temporal coherence {coherence_p2:.2f}.")
    cand_mask2 = selectPixels(
        path=output_path, selection_method="temp_coh",
        thrsh=coherence_p2,
        grid_size=None, bool_plot=True,
        logger=logger
    )  # first-order points are included in second-order points
    logger.debug(f"{cand_mask2.sum()} second-order points selected based on temporal coherence {coherence_p2:.2f}.")

    if config.phase_linking.phase_linking:
        # read PL results
        pl_file = join(config.phase_linking.path_inverted, "phase_series.h5")
        logger.debug(f"Read phase linking results from file {pl_file}.")
        pl_coh = readfile.read(pl_file,
                               datasetName='temporalCoherence')[0]
        pl_coh = pl_coh[1, :, :]
        siblings = readfile.read(pl_file,
                                 datasetName='shp')[0]

        if config.phase_linking.use_ps:
            logger.debug(f"Read phase linking PS mask from file {config.phase_linking.path_mask_file_ps}.")
            mask_ps = readfile.read(config.phase_linking.path_mask_file_ps,
                                    datasetName='mask')[0].astype(np.bool_)
            cand_mask_pl = (pl_coh > coherence_p2) | mask_ps
        else:
            cand_mask_pl = (pl_coh > coherence_p2)
            # remove ps candidates, because the ps detection strategy in miaplpy seems to be biased.
            cand_mask_pl[siblings <= config.phase_linking.num_siblings] = False

        if config.phase_linking.spatial_mask_file_pl is not None:
            path_mask_pl_aoi = join(config.phase_linking.spatial_mask_file_pl)
            logger.debug(f"Load mask for area of interest from file {path_mask_pl_aoi}.")
            mask_pl_aoi = readfile.read(path_mask_pl_aoi, datasetName='mask')[0].astype(np.bool_)

            fig = plt.figure(figsize=[15, 5])
            ax = fig.add_subplot()
            ax.imshow(mask_pl_aoi, cmap=plt.cm.get_cmap("gray"), alpha=0.5, zorder=10, vmin=0, vmax=1)
            bmap_obj.plot(ax=ax, logger=logger)
            coord_xy = np.array(np.where(cand_mask_pl)).transpose()
            val = np.ones_like(cand_mask_pl)
            sc = ax.scatter(coord_xy[:, 1], coord_xy[:, 0], c=val[cand_mask_pl], s=0.5,
                            cmap=plt.get_cmap("autumn_r"),
                            vmin=1, vmax=2)  # set min, max to ensure that points are yellow
            cbar = plt.colorbar(sc, pad=0.03, shrink=0.5)
            cbar.ax.set_visible(False)  # make size of axis consistent with all others
            plt.tight_layout()
            plt.title("Mask for phase linking points")
            fig.savefig(join(output_path, "pic", "step_3_mask_coh{}_phase_linking.png".format(coh_value)), dpi=300)
            plt.close(fig)

            # mask points after plotting, so that available coherent points are visible in figure
            cand_mask_pl[~mask_pl_aoi] = False

        # combine phase linking coherence with TPC cand_mask2
        initial_num_cand_mask2 = cand_mask2.sum()
        cand_mask2 = cand_mask2 | cand_mask_pl
        logger.debug(f"{cand_mask2.sum()-initial_num_cand_mask2} additional pixels selected by phase linking.")

    mask_valid_area = ut.detectValidAreas(bmap_obj=bmap_obj, logger=logger)

    if config.filtering.spatial_mask_file_p2 is not None:
        path_mask_aoi = join(config.filtering.spatial_mask_file_p2)
        logger.info(f"Load mask for area of interest from file {path_mask_aoi}.")
        mask_aoi = readfile.read(path_mask_aoi, datasetName='mask')[0].astype(np.bool_)
        mask_valid_area &= mask_aoi
        # todo: add unstable points from p1 for densification
    else:
        logger.info(msg="No mask for area of interest given.")

    cand_mask2[p1_mask] = True  # add all selected 1.order points to avoid spatial gaps in 2D unwrapping
    # cand_mask2[cand_mask_sparse] = True  # add only stable points from 1.order points

    cand_mask2 &= mask_valid_area
    logger.info(f"Final number of selected second-order points: {cand_mask2.sum()}.")

    fig = plt.figure(figsize=[15, 5])
    ax = fig.add_subplot()
    ax.imshow(mask_valid_area, cmap=plt.cm.get_cmap("gray"), alpha=0.5, zorder=10, vmin=0, vmax=1)
    bmap_obj.plot(ax=ax, logger=logger)
    coord_xy = np.array(np.where(cand_mask2)).transpose()
    val = np.ones_like(cand_mask2)
    sc = ax.scatter(coord_xy[:, 1], coord_xy[:, 0], c=val[cand_mask2], s=0.5, cmap=plt.get_cmap("autumn_r"),
                    vmin=1, vmax=2)  # set min, max to ensure that points are yellow
    cbar = plt.colorbar(sc, pad=0.03, shrink=0.5)
    cbar.ax.set_visible(False)  # make size of axis consistent with all others
    plt.tight_layout()
    plt.title("Mask for dense point set")
    fig.savefig(join(output_path, "pic", "step_3_mask_coh{}.png".format(coh_value)), dpi=300)
    plt.close(fig)

    p2_file = join(output_path, f"coh{coh_value}_ifg_wr.h5")
    logger.debug(f"Prepare second-order points object to save to file {p2_file}.")
    point2_obj = Points(file_path=p2_file, logger=logger)
    coord_xy = np.array(np.where(cand_mask2)).transpose()
    point_id2 = point_id_img[cand_mask2]
    point2_obj.prepare(
        point_id=point_id2,
        coord_xy=coord_xy,
        path_inputs=config.data_directories.path_inputs
    )

    ifg_stack_file = join(output_path, "ifg_stack.h5")
    logger.debug(f"Read interferogram phase for selected second-order points from file {ifg_stack_file}.")
    ifg_stack_obj = BaseStack(file=ifg_stack_file, logger=logger)
    point2_obj.phase = ut.readPhasePatchwise(stack_obj=ifg_stack_obj, dataset_name="ifgs",
                                             num_patches=config.processing.num_patches, cand_mask=cand_mask2,
                                             point_id_img=point_id_img, logger=logger)

    if config.phase_linking.phase_linking:
        pl_inverted_file = join(config.phase_linking.path_inverted, "phase_series.h5")
        logger.debug(f"Read interferogram phase from MiaplPy results from file {pl_inverted_file}.")
        phase_linking_obj = BaseStack(file=pl_inverted_file,
                                      logger=logger)

        pl_phase = ut.readPhasePatchwise(
            stack_obj=phase_linking_obj, dataset_name="phase",
            num_patches=config.processing.num_patches,
            cand_mask=cand_mask2,
            point_id_img=point_id_img, logger=logger
        )

        # subset to time span
        slc_stack_obj = slcStack(join(config.data_directories.path_inputs, "slcStack.h5"))
        slc_stack_obj.open()
        time_mask = createTimeMaskFromDates(
            start_date=config.preparation.start_date,
            stop_date=config.preparation.stop_date,
            date_list=slc_stack_obj.dateList,
            logger=logger
        )[0]
        pl_phase = pl_phase[:, time_mask]

        pl_ifgs = np.zeros((point2_obj.num_points, point2_obj.ifg_net_obj.num_ifgs), dtype=np.float32)

        c = 0
        for i, j in np.asarray(point1_obj.ifg_net_obj.ifg_list):
            pl_ifgs[:, c] = np.angle(np.exp(1j * pl_phase[:, i]) * np.conjugate(np.exp(1j * pl_phase[:, j])))
            c += 1

        # change only phase for good phase linking pixels and keep original phase for good tpc pixels
        mask_pl = cand_mask_pl[cand_mask2]
        point2_obj.phase[mask_pl] = pl_ifgs[mask_pl]

    logger.debug(f"Write second-order points to file {point2_obj.file_path}.")
    point2_obj.writeToFile()
    del ifg_stack_obj

    p1_aps_file = join(output_path, "p1_aps.h5")
    logger.debug(f"Read APS for first order point from file {p1_aps_file}.")
    aps1_obj = Points(file_path=p1_aps_file, logger=logger)
    aps1_obj.open(path_inputs=config.data_directories.path_inputs)

    aps_params_file = join(output_path, "aps_parameters.h5")
    logger.debug(f"Read APS parameters from file {aps_params_file}.")
    aps_params_obj = ApsParameters(file_path=aps_params_file, logger=logger)
    aps_params_obj.open()

    if config.filtering.skip_filtering:
        logger.info("Skip APS filtering")
        # TODO: check if this is correct ####
        # original code: aps2_phase = np.zeros_like(point2_obj.phase)
        aps2_phase = np.zeros((point2_obj.num_points, aps_params_obj.phase.shape[1]))
    else:
        logger.debug("Apply APS filtering to second-order points.")
        if config.filtering.interpolation_method == "kriging":
            aps2_phase = applySpatialFilteringToP2(coord_utm1=aps1_obj.coord_utm,
                                                   residuals=aps_params_obj.phase,
                                                   coord_utm2=point2_obj.coord_utm,
                                                   model_name=aps_params_obj.model_name,
                                                   model_params=aps_params_obj.model_params,
                                                   logger=logger)
        else:
            aps2_phase = applySimpleInterpolationToP2(residuals=aps_params_obj.phase,
                                                      coord_utm1=aps1_obj.coord_utm,
                                                      coord_utm2=point2_obj.coord_utm,
                                                      logger=logger,
                                                      interp_method=config.filtering.interpolation_method)

    p2_aps_file = join(output_path, f"coh{coh_value}_aps.h5")
    logger.info(f"Write APS for second-order points to file {p2_aps_file}.")
    aps2_obj = Points(file_path=p2_aps_file, logger=logger)
    aps2_obj.open(
        other_file_path=join(output_path, f"coh{coh_value}_ifg_wr.h5"),
        path_inputs=config.data_directories.path_inputs
    )
    aps2_obj.phase = aps2_phase
    aps2_obj.writeToFile()

    logger.debug("Finished selecting second-order points.")
    return point2_obj, aps2_obj
