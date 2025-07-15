#!/usr/bin/env python

# SARvey - A multitemporal InSAR time series tool for the derivation of displacements.
#
# Copyright (C) 2021-2025 Andreas Piter (IPI Hannover, piter@ipi.uni-hannover.de)
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

"""Processing module for SARvey."""
from os.path import join
import matplotlib.pyplot as plt
import numpy as np
from logging import Logger
import cmcrameri as cmc

from miaplpy.objects.slcStack import slcStack
from mintpy.utils import readfile
from mintpy.utils.plot import auto_flip_direction

from sarvey import viewer
from sarvey.densification import densifyNetwork
from sarvey.filtering import estimateAtmosphericPhaseScreen, simpleInterpolation
from sarvey.ifg_network import (DelaunayNetwork, SmallBaselineYearlyNetwork, SmallTemporalBaselinesNetwork,
                                SmallBaselineNetwork, StarNetwork)
from sarvey.objects import Network, Points, AmplitudeImage, CoordinatesMap, NetworkParameter, BaseStack
from sarvey.unwrapping import spatialParameterIntegration, \
    parameterBasedNoisyPointRemoval, temporalUnwrapping, spatialUnwrapping, removeGrossOutliers
from sarvey.preparation import createArcsBetweenPoints, selectPixels, createTimeMaskFromDates
import sarvey.utils as ut
from sarvey.coherence import computeIfgsAndTemporalCoherence
from sarvey.triangulation import PointNetworkTriangulation
from sarvey.config import Config


class Processing:
    """Processing."""

    def __init__(self, path: str, config: Config, logger: Logger):
        """Init."""
        self.path = path
        self.config = config
        self.logger = logger

    def runPreparation(self):
        """RunPreparation."""
        log = self.logger

        msg = "#" * 10
        msg += " PREPARE PROCESSING: LOAD INPUT "
        msg += "#" * 10
        log.info(msg=msg)

        # load slc data
        slc_stack_obj = slcStack(join(self.config.general.input_path, "slcStack.h5"))
        slc_stack_obj.open()

        if "ORBIT_DIRECTION" in slc_stack_obj.metadata:
            log.info(msg="Orbit direction: {}".format(slc_stack_obj.metadata["ORBIT_DIRECTION"]))
        else:
            log.warning(msg="No orbit direction found in metadata. Add 'ORBIT_DIRECTION' to metadata of 'slcStack.h5'"
                            "and run again!")
            raise AttributeError("No orbit direction found in metadata.")

        time_mask, num_slc, date_list = createTimeMaskFromDates(
            start_date=self.config.preparation.start_date,
            stop_date=self.config.preparation.end_date,
            date_list=slc_stack_obj.dateList,
            logger=log
        )
        log.info(msg=f"Start date: {date_list[0]}")
        log.info(msg=f"Stop date: {date_list[-1]}")
        log.info(msg=f"Number of SLC: {num_slc}")

        msg = "#" * 10
        msg += " DESIGN IFG NETWORK "
        msg += "#" * 10
        log.info(msg=msg)

        ifg_net_obj = None
        if self.config.preparation.ifg_network_type == "star":
            ifg_net_obj = StarNetwork()
            ifg_net_obj.configure(
                pbase=slc_stack_obj.pbase[time_mask],
                tbase=slc_stack_obj.tbase[time_mask],
                ref_idx=int(np.floor(num_slc/2)),
                dates=date_list
            )
            log.info(msg="Star ifg network")
        elif self.config.preparation.ifg_network_type == "sb":
            ifg_net_obj = SmallBaselineNetwork()
            ifg_net_obj.configure(pbase=slc_stack_obj.pbase[time_mask],
                                  tbase=slc_stack_obj.tbase[time_mask],
                                  num_link=self.config.preparation.num_ifgs,
                                  max_tbase=self.config.preparation.max_tbase,
                                  dates=date_list)
            log.info(msg="Small baseline network")
        elif self.config.preparation.ifg_network_type == "stb":
            ifg_net_obj = SmallTemporalBaselinesNetwork()
            ifg_net_obj.configure(
                pbase=slc_stack_obj.pbase[time_mask],
                tbase=slc_stack_obj.tbase[time_mask],
                num_link=self.config.preparation.num_ifgs,
                dates=date_list
            )
            log.info(msg="Small temporal baseline network")
        elif self.config.preparation.ifg_network_type == "stb_year":
            ifg_net_obj = SmallBaselineYearlyNetwork()
            ifg_net_obj.configure(
                pbase=slc_stack_obj.pbase[time_mask],
                tbase=slc_stack_obj.tbase[time_mask],
                num_link=self.config.preparation.num_ifgs,
                dates=date_list
            )
            log.info(msg="Small temporal baseline and yearly ifg network")
        elif self.config.preparation.ifg_network_type == "delaunay":
            ifg_net_obj = DelaunayNetwork()
            ifg_net_obj.configure(
                pbase=slc_stack_obj.pbase[time_mask],
                tbase=slc_stack_obj.tbase[time_mask],
                dates=date_list
            )
            log.info(msg="Delaunay ifg network")

        ifg_net_obj.writeToFile(path=join(self.path, "ifg_network.h5"), logger=log)
        log.info(msg=f"temporal baselines: {np.unique(np.round(np.abs(ifg_net_obj.tbase_ifg) * 365.25).astype(int))}")

        fig = ifg_net_obj.plot()
        fig.savefig(join(self.path, "pic", "step_0_interferogram_network.png"), dpi=300)
        plt.close(fig)
        # at this point just created folder pic and ifg_network.h5
        msg = "#" * 10
        msg += f" GENERATE STACK OF {ifg_net_obj.num_ifgs} INTERFEROGRAMS & ESTIMATE TEMPORAL COHERENCE "
        msg += "#" * 10
        log.info(msg=msg)

        box_list, num_patches = ut.preparePatches(num_patches=self.config.general.num_patches,
                                                  width=slc_stack_obj.width,
                                                  length=slc_stack_obj.length,
                                                  logger=log)

        # create placeholder in result file for datasets which are stored patch-wise
        dshape = (slc_stack_obj.length, slc_stack_obj.width, ifg_net_obj.num_ifgs)
        ifg_stack_obj = BaseStack(file=join(self.path, "ifg_stack.h5"), logger=log)
        ifg_stack_obj.prepareDataset(dataset_name="ifgs", dshape=dshape, dtype=np.csingle,
                                     metadata=slc_stack_obj.metadata, mode='w', chunks=(30, 30, ifg_net_obj.num_ifgs))

        # create placeholder in result file for datasets which are stored patch-wise
        temp_coh_obj = BaseStack(file=join(self.path, "temporal_coherence.h5"), logger=log)
        dshape = (slc_stack_obj.length, slc_stack_obj.width)
        temp_coh_obj.prepareDataset(dataset_name="temp_coh", metadata=slc_stack_obj.metadata,
                                    dshape=dshape, dtype=np.float32, mode="w", chunks=True)

        mean_amp_img = computeIfgsAndTemporalCoherence(
            path_temp_coh=join(self.path, "temporal_coherence.h5"),
            path_ifgs=join(self.path, "ifg_stack.h5"),
            path_slc=join(self.config.general.input_path, "slcStack.h5"),
            ifg_array=np.array(ifg_net_obj.ifg_list),
            time_mask=time_mask,
            wdw_size=self.config.preparation.filter_window_size,
            num_boxes=num_patches,
            box_list=box_list,
            num_cores=self.config.general.num_cores,
            logger=log)

        # store auxilliary datasets for faster access during processing
        coord_map_obj = CoordinatesMap(file_path=join(self.path, "coordinates_map.h5"), logger=self.logger)
        coord_map_obj.prepare(input_path=join(self.config.general.input_path, "geometryRadar.h5"))
        del coord_map_obj

        bmap_obj = AmplitudeImage(file_path=join(self.path, "background_map.h5"))
        bmap_obj.prepare(slc_stack_obj=slc_stack_obj, img=mean_amp_img, logger=self.logger)
        ax = bmap_obj.plot(logger=self.logger)
        img = ax.get_images()[0]
        cbar = plt.colorbar(img, pad=0.03, shrink=0.5)
        cbar.ax.set_visible(False)
        plt.tight_layout()
        plt.gcf().savefig(join(self.path, "pic", "step_0_amplitude_image.png"), dpi=300)
        plt.close(plt.gcf())
        del bmap_obj
        del mean_amp_img

        temp_coh = temp_coh_obj.read(dataset_name="temp_coh")

        fig = plt.figure(figsize=(15, 5))
        ax = fig.add_subplot()
        im = ax.imshow(temp_coh, cmap=cmc.cm.cmaps["grayC"], vmin=0, vmax=1)
        auto_flip_direction(slc_stack_obj.metadata, ax=ax, print_msg=True)
        ax.set_xlabel("Range")
        ax.set_ylabel("Azimuth")
        plt.colorbar(im, pad=0.03, shrink=0.5)
        plt.title("Temporal coherence")
        plt.tight_layout()
        fig.savefig(join(self.path, "pic", "step_0_temporal_phase_coherence.png"), dpi=300)
        plt.close(fig)

    def runConsistencyCheck(self):
        """RunConsistencyCheck."""
        # 0) select candidates for first order points
        ifg_stack_obj = BaseStack(file=join(self.path, "ifg_stack.h5"), logger=self.logger)
        length, width, num_ifgs = ifg_stack_obj.getShape(dataset_name="ifgs")

        cand_mask1 = selectPixels(
            path=self.path, selection_method="temp_coh", thrsh=self.config.consistency_check.coherence_p1,
            grid_size=self.config.consistency_check.grid_size, bool_plot=True, logger=self.logger
        )

        bmap_obj = AmplitudeImage(file_path=join(self.path, "background_map.h5"))
        mask_valid_area = ut.detectValidAreas(bmap_obj=bmap_obj, logger=self.logger)

        if self.config.consistency_check.mask_p1_file is not None:
            path_mask_aoi = join(self.config.consistency_check.mask_p1_file)
            self.logger.info(msg="load mask for area of interest from: {}.".format(path_mask_aoi))
            mask_aoi = readfile.read(path_mask_aoi, datasetName='mask')[0].astype(np.bool_)
            mask_valid_area &= mask_aoi
        else:
            self.logger.info(msg="No mask for area of interest given.")

        cand_mask1 &= mask_valid_area

        fig = plt.figure(figsize=(15, 5))
        ax = fig.add_subplot()
        ax.imshow(mask_valid_area, cmap=cmc.cm.cmaps["grayC"], alpha=0.5, zorder=10, vmin=0, vmax=1)
        bmap_obj.plot(ax=ax, logger=self.logger)
        coord_xy = np.array(np.where(cand_mask1)).transpose()
        val = np.ones_like(cand_mask1)
        sc = ax.scatter(coord_xy[:, 1], coord_xy[:, 0], c=val[cand_mask1], s=0.5, cmap=cmc.cm.cmaps["lajolla_r"],
                        vmin=1, vmax=2)  # set min, max to ensure that points are yellow
        cbar = plt.colorbar(sc, pad=0.03, shrink=0.5)
        cbar.ax.set_visible(False)  # make size of axis consistent with all others
        plt.tight_layout()
        plt.title("Mask for first order point set")
        fig.savefig(join(self.path, "pic", "step_1_mask_p1.png"), dpi=300)
        plt.close(fig)

        if cand_mask1[cand_mask1].shape[0] == 0:
            self.logger.error("No points selected for first-order points. Modify the coherence threshold.")
            raise ValueError

        # create unique point_id throughout the image to make it possible to mix first-order and second-order points
        # in the densification step. point_id is ordered so that it fits to anydata[mask].ravel() when loading the data.
        point_id_img = np.arange(0, length * width).reshape((length, width))

        point_obj = Points(file_path=join(self.path, "p1_ifg_wr.h5"), logger=self.logger)
        point_id1 = point_id_img[cand_mask1]

        point_obj.prepare(
            point_id=point_id1,
            coord_xy=coord_xy,
            input_path=self.config.general.input_path
        )

        point_obj.phase = ut.readPhasePatchwise(stack_obj=ifg_stack_obj, dataset_name="ifgs",
                                                num_patches=self.config.general.num_patches, cand_mask=cand_mask1,
                                                point_id_img=point_id_img, logger=self.logger)

        point_obj.writeToFile()
        del ifg_stack_obj, cand_mask1

        # 1) create spatial network
        arcs = createArcsBetweenPoints(point_obj=point_obj,
                                       knn=self.config.consistency_check.num_nearest_neighbours,
                                       max_arc_length=self.config.consistency_check.max_arc_length,
                                       logger=self.logger)
        net_obj = Network(file_path=join(self.path, "point_network.h5"), logger=self.logger)
        net_obj.computeArcObservations(
            point_obj=point_obj,
            arcs=arcs
        )
        net_obj.writeToFile()
        net_obj.open(input_path=self.config.general.input_path)  # to retrieve external data

        demerr, vel, gamma = temporalUnwrapping(ifg_net_obj=point_obj.ifg_net_obj,
                                                net_obj=net_obj,
                                                wavelength=point_obj.wavelength,
                                                velocity_bound=self.config.consistency_check.velocity_bound,
                                                demerr_bound=self.config.consistency_check.dem_error_bound,
                                                num_samples=self.config.consistency_check.num_optimization_samples,
                                                num_cores=self.config.general.num_cores,
                                                logger=self.logger)

        net_par_obj = NetworkParameter(file_path=join(self.path, "point_network_parameter.h5"),
                                       logger=self.logger)
        net_par_obj.prepare(
            net_obj=net_obj,
            demerr=demerr,
            vel=vel,
            gamma=gamma
        )
        net_par_obj.writeToFile()

        # 3) spatial unwrapping of the arc network and removal of outliers (arcs and points)
        bmap_obj = AmplitudeImage(file_path=join(self.path, "background_map.h5"))
        thrsh_visualisation = 0.7

        try:
            ax = bmap_obj.plot(logger=self.logger)
            arc_mask = net_par_obj.gamma.reshape(-1) <= thrsh_visualisation
            ax, cbar = viewer.plotColoredPointNetwork(x=point_obj.coord_xy[:, 1], y=point_obj.coord_xy[:, 0],
                                                      arcs=net_par_obj.arcs[arc_mask, :],
                                                      val=net_par_obj.gamma[arc_mask],
                                                      ax=ax, linewidth=1, cmap="lajolla", clim=(0, 1))
            ax.set_title("Coherence from temporal unwrapping\n"
                         r"(only arcs with $\gamma \leq$ {} "
                         "shown)\nBefore outlier removal".format(thrsh_visualisation))
            fig = ax.get_figure()
            plt.tight_layout()
            fig.savefig(join(self.path, "pic", "step_1_arc_coherence.png"), dpi=300)
        except BaseException as e:
            self.logger.exception(msg="NOT POSSIBLE TO PLOT SPATIAL NETWORK OF POINTS. {}".format(e))

        net_par_obj, point_id, coord_xy, design_mat = removeGrossOutliers(
            net_obj=net_par_obj,
            point_id=point_obj.point_id,
            coord_xy=point_obj.coord_xy,
            min_num_arc=self.config.consistency_check.min_num_arc,
            quality_thrsh=self.config.consistency_check.arc_unwrapping_coherence,
            logger=self.logger
        )

        try:
            ax = bmap_obj.plot(logger=self.logger)
            arc_mask = net_par_obj.gamma.reshape(-1) <= thrsh_visualisation
            ax, cbar = viewer.plotColoredPointNetwork(x=coord_xy[:, 1], y=coord_xy[:, 0],
                                                      arcs=net_par_obj.arcs[arc_mask, :],
                                                      val=net_par_obj.gamma[arc_mask],
                                                      ax=ax, linewidth=1, cmap="lajolla", clim=(0, 1))
            ax.set_title("Coherence from temporal unwrapping\n"
                         r"(only arcs with $\gamma \leq$ {} "
                         "shown)\nAfter outlier removal".format(thrsh_visualisation))
            fig = ax.get_figure()
            plt.tight_layout()
            fig.savefig(join(self.path, "pic", "step_1_arc_coherence_reduced.png"), dpi=300)
        except BaseException as e:
            self.logger.exception(msg="NOT POSSIBLE TO PLOT SPATIAL NETWORK OF POINTS. {}".format(e))

        spatial_ref_id, point_id, net_par_obj = parameterBasedNoisyPointRemoval(
            net_par_obj=net_par_obj,
            point_id=point_id,
            coord_xy=coord_xy,
            design_mat=design_mat,
            bmap_obj=bmap_obj,
            bool_plot=True,
            logger=self.logger
        )

        net_par_obj.writeToFile()  # arcs were removed. obj still needed in next step.
        point_obj.removePoints(keep_id=point_id, input_path=self.config.general.input_path)
        point_obj.writeToFile()

    def runUnwrappingTimeAndSpace(self):
        """RunTemporalAndSpatialUnwrapping."""
        net_par_obj = NetworkParameter(file_path=join(self.path, "point_network_parameter.h5"),
                                       logger=self.logger)
        net_par_obj.open(input_path=self.config.general.input_path)

        point_obj = Points(file_path=join(self.path, "p1_ifg_unw.h5"), logger=self.logger)
        point_obj.open(
            other_file_path=join(self.path, "p1_ifg_wr.h5"),
            input_path=self.config.general.input_path
        )

        # reference point can be set arbitrarily, because outliers are removed.
        spatial_ref_idx = 0

        bmap_obj = AmplitudeImage(file_path=join(self.path, "background_map.h5"))

        self.logger.info(msg="Integrate parameters from arcs to points.")
        self.logger.info(msg="Integrate DEM correction.")
        demerr = spatialParameterIntegration(val_arcs=net_par_obj.demerr,
                                             arcs=net_par_obj.arcs,
                                             coord_xy=point_obj.coord_xy,
                                             weights=net_par_obj.gamma,
                                             spatial_ref_idx=spatial_ref_idx, logger=self.logger)

        # demerr = spatialParameterIntegrationIterative(val_arcs=net_par_obj.demerr, all_arcs=net_par_obj.arcs,
        #                                               coord_xy=point_obj.coord_xy, all_weights=net_par_obj.gamma,
        #                                               spatial_ref_idx=spatial_ref_idx,
        #                                               res_tol=5.0,
        #                                               max_rm_fraction=0.001)
        fig = viewer.plotScatter(value=-demerr, coord=point_obj.coord_xy,
                                 ttl="Parameter integration: DEM correction in [m]",
                                 bmap_obj=bmap_obj, s=3.5, cmap="vanimo", symmetric=True,
                                 logger=self.logger)[0]
        fig.savefig(join(self.path, "pic", "step_2_estimation_dem_correction.png"), dpi=300)
        plt.close(fig)

        self.logger.info(msg="Integrate mean velocity.")
        vel = spatialParameterIntegration(val_arcs=net_par_obj.vel,
                                          arcs=net_par_obj.arcs,
                                          coord_xy=point_obj.coord_xy,
                                          weights=net_par_obj.gamma,
                                          spatial_ref_idx=spatial_ref_idx, logger=self.logger)

        # vel = spatialParameterIntegrationIterative(val_arcs=net_par_obj.vel, all_arcs=net_par_obj.arcs,
        #                                            coord_xy=point_obj.coord_xy,
        #                                            all_weights=net_par_obj.gamma,
        #                                            spatial_ref_idx=spatial_ref_idx,
        #                                            res_tol=1.0,
        #                                            max_rm_fraction=0.001)
        fig = viewer.plotScatter(value=-vel, coord=point_obj.coord_xy,
                                 ttl="Parameter integration: mean velocity in [m / year]",
                                 bmap_obj=bmap_obj, s=3.5, cmap="roma", symmetric=True,
                                 logger=self.logger)[0]
        fig.savefig(join(self.path, "pic", "step_2_estimation_velocity.png"), dpi=300)
        plt.close(fig)

        self.logger.info(msg="Remove phase contributions from mean velocity"
                             " and DEM correction from wrapped phase of points.")
        pred_phase_demerr, pred_phase_vel = ut.predictPhase(
            obj=point_obj, vel=vel, demerr=demerr,
            ifg_space=True, logger=self.logger
        )
        pred_phase = pred_phase_demerr + pred_phase_vel

        wr_phase = point_obj.phase
        wr_res_phase = np.angle(np.exp(1j * wr_phase) * np.conjugate(np.exp(1j * pred_phase)))

        if self.config.unwrapping.use_arcs_from_temporal_unwrapping:
            arcs = net_par_obj.arcs  # use this to avoid unreliable connections. Takes a bit longer.
        else:
            triang_obj = PointNetworkTriangulation(coord_xy=point_obj.coord_xy, coord_map_xy=point_obj.coord_map,
                                                   logger=self.logger)
            triang_obj.triangulateGlobal()
            arcs = triang_obj.getArcsFromAdjMat()

        unw_res_phase = spatialUnwrapping(num_ifgs=point_obj.ifg_net_obj.num_ifgs,
                                          num_points=point_obj.num_points,
                                          phase=wr_res_phase,
                                          method=self.config.general.spatial_unwrapping_method,
                                          edges=arcs,
                                          num_cores=self.config.general.num_cores, logger=self.logger)

        # use same reference point for spatial integration and Puma unwrapping before recombining phases
        unw_res_phase = unw_res_phase - unw_res_phase[spatial_ref_idx, :]

        self.logger.info(msg="Add phase contributions from mean velocity and DEM correction back to "
                             "spatially unwrapped residual phase.")
        unw_phase = unw_res_phase + pred_phase
        # unw_phase = unw_res_phase  # debug: don't add phase back.

        # adjust reference to peak of histogram
        point_obj.phase = unw_phase
        vel = ut.estimateParameters(obj=point_obj, ifg_space=True)[0]
        point_obj.phase = ut.setReferenceToPeakOfHistogram(phase=unw_phase, vel=vel, num_bins=300)

        point_obj.writeToFile()

        phase_ts = ut.invertIfgNetwork(
            phase=unw_phase,
            num_points=point_obj.num_points,
            ifg_net_obj=point_obj.ifg_net_obj,
            logger=self.logger
        )
        point_obj = Points(file_path=join(self.path, "p1_ts.h5"), logger=self.logger)
        point_obj.open(
            other_file_path=join(self.path, "p1_ifg_unw.h5"),
            input_path=self.config.general.input_path
        )
        point_obj.phase = phase_ts
        point_obj.writeToFile()

    def runUnwrappingSpace(self):
        """RunSpatialUnwrapping."""
        point_obj = Points(file_path=join(self.path, "p1_ifg_unw.h5"), logger=self.logger)
        point_obj.open(
            other_file_path=join(self.path, "p1_ifg_wr.h5"),
            input_path=self.config.general.input_path
        )

        if self.config.unwrapping.use_arcs_from_temporal_unwrapping:
            net_par_obj = NetworkParameter(file_path=join(self.path, "point_network_parameter.h5"),
                                           logger=self.logger)
            net_par_obj.open(input_path=self.config.general.input_path)
            arcs = net_par_obj.arcs  # use this to avoid unreliable connections. Takes a bit longer.
        else:
            # re-triangulate with delaunay to make PUMA faster
            triang_obj = PointNetworkTriangulation(coord_xy=point_obj.coord_xy, coord_map_xy=point_obj.coord_map,
                                                   logger=self.logger)
            triang_obj.triangulateGlobal()
            arcs = triang_obj.getArcsFromAdjMat()

        bmap_obj = AmplitudeImage(file_path=join(self.path, "background_map.h5"))
        ax = bmap_obj.plot(logger=self.logger)
        ax, cbar = viewer.plotColoredPointNetwork(x=point_obj.coord_xy[:, 1],
                                                  y=point_obj.coord_xy[:, 0],
                                                  arcs=arcs,
                                                  val=np.zeros(arcs.shape[0], dtype=np.float32),
                                                  ax=ax, linewidth=0.5, cmap="lajolla", clim=(0, 1))
        cbar.ax.set_visible(False)
        ax.set_xlabel("Range")
        ax.set_ylabel("Azimuth")
        ax.set_title("Unwrapping Network")
        plt.tight_layout()
        plt.gcf().savefig(join(self.path, "pic", "step_2_unwrapping_network_p1.png"), dpi=300)
        plt.close(plt.gcf())

        unw_phase = spatialUnwrapping(num_ifgs=point_obj.ifg_net_obj.num_ifgs,
                                      num_points=point_obj.num_points,
                                      phase=point_obj.phase,
                                      method=self.config.general.spatial_unwrapping_method,
                                      edges=arcs,
                                      num_cores=self.config.general.num_cores, logger=self.logger)

        # adjust reference to peak of histogram
        point_obj.phase = unw_phase
        vel = ut.estimateParameters(obj=point_obj, ifg_space=True)[0]
        point_obj.phase = ut.setReferenceToPeakOfHistogram(phase=unw_phase, vel=vel, num_bins=300)

        point_obj.writeToFile()
        del point_obj

        point_obj = Points(file_path=join(self.path, "p1_ts.h5"), logger=self.logger)
        point_obj.open(
            other_file_path=join(self.path, "p1_ifg_wr.h5"),
            input_path=self.config.general.input_path
        )

        # for sbas the ifg network needs to be inverted to get the phase time series
        phase_ts = ut.invertIfgNetwork(phase=unw_phase, num_points=point_obj.num_points,
                                       ifg_net_obj=point_obj.ifg_net_obj,
                                       logger=self.logger)

        point_obj.phase = phase_ts
        point_obj.writeToFile()

    def runFiltering(self):
        """RunFiltering."""
        coh_value = int(self.config.filtering.coherence_p2 * 100)

        # create output file which contains filtered phase time series
        point1_obj = Points(file_path=join(self.path, "p1_ts_filt.h5"), logger=self.logger)
        point1_obj.open(
            other_file_path=join(self.path, "p1_ts.h5"),
            input_path=self.config.general.input_path
        )
        p1_mask = point1_obj.createMask()  # used later for selecting psCand2 when a spatial mask AOI is given.

        # select only pixels which have low phase noise and are well distributed
        mask = point1_obj.createMask()

        bmap_obj = AmplitudeImage(file_path=join(self.path, "background_map.h5"))

        # temporal auto-correlation
        auto_corr_img = np.zeros_like(mask, np.float64)

        vel, demerr, _, _, _, residuals = ut.estimateParameters(obj=point1_obj, ifg_space=False)

        if self.config.filtering.use_moving_points:
            auto_corr = ut.temporalAutoCorrelation(residuals=residuals, lag=1).reshape(-1)
        else:
            # remove DEM error, but not velocity before estimating the temporal autocorrelation
            pred_phase_demerr = ut.predictPhase(
                obj=point1_obj, vel=vel, demerr=demerr, ifg_space=False, logger=self.logger)[0]
            phase_wo_demerr = point1_obj.phase - pred_phase_demerr
            auto_corr = ut.temporalAutoCorrelation(residuals=phase_wo_demerr, lag=1).reshape(-1)

        auto_corr_img[mask] = auto_corr
        auto_corr_img[~mask] = np.inf

        fig = viewer.plotScatter(value=auto_corr, coord=point1_obj.coord_xy, bmap_obj=bmap_obj,
                                 ttl="Temporal autocorrelation", unit="[ ]", s=3.5, cmap="lajolla",
                                 vmin=0, vmax=1, logger=self.logger)[0]
        fig.savefig(join(self.path, "pic", "step_3_temporal_autocorrelation.png"), dpi=300)
        plt.close(fig)

        # create grid
        coord_map_obj = CoordinatesMap(file_path=join(self.path, "coordinates_map.h5"), logger=self.logger)
        coord_map_obj.open()

        # remove points based on threshold
        mask_thrsh = auto_corr_img <= self.config.filtering.max_temporal_autocorrelation
        auto_corr_img[~mask_thrsh] = np.inf

        box_list, num_box = ut.createSpatialGrid(coord_map_img=coord_map_obj.coord_map, length=point1_obj.length,
                                                 width=point1_obj.width,
                                                 grid_size=self.config.filtering.grid_size,
                                                 logger=self.logger)

        cand_mask_sparse = ut.selectBestPointsInGrid(box_list=box_list, quality=auto_corr_img, sel_min=True)

        num_p1_points_for_filtering = cand_mask_sparse[cand_mask_sparse].shape[0]
        if num_p1_points_for_filtering < 10:
            self.logger.warning(msg=f"Only {num_p1_points_for_filtering} points for APS filtering selected. Filtering "
                                    f"results are probably not reliable. You can e.g. increase 'max_auto_corr' or try "
                                    f"to increase the number of first-order points during step 1 and 2.")
            if num_p1_points_for_filtering == 0:
                self.logger.error("No points selected for APS filtering.")
                raise ValueError

        point_id_img = np.arange(0, point1_obj.length * point1_obj.width).reshape(
            (point1_obj.length, point1_obj.width))
        keep_id = point_id_img[np.where(cand_mask_sparse)]
        point1_obj.removePoints(keep_id=keep_id, input_path=self.config.general.input_path)
        point1_obj.writeToFile()  # to be able to load aps1 from this file having the same set of points

        fig, ax = viewer.plotScatter(value=auto_corr_img[cand_mask_sparse], coord=point1_obj.coord_xy,
                                     bmap_obj=bmap_obj, ttl="Selected pixels for APS estimation",
                                     unit="Auto-correlation\n[ ]", s=5, cmap="lajolla", vmin=0, vmax=1,
                                     logger=self.logger)[:2]
        viewer.plotGridFromBoxList(box_list=box_list, ax=ax, edgecolor="k", linewidth=0.2)
        fig.savefig(join(self.path, "pic", "step_3_stable_points.png"), dpi=300)
        plt.close(fig)

        if self.config.filtering.use_moving_points:
            # recompute the residuals, because now there are fewer points in the obj
            phase_for_aps_filtering = ut.estimateParameters(obj=point1_obj, ifg_space=False)[-1]
        else:
            phase_for_aps_filtering = point1_obj.phase

        # create output which contains only the atmospheric phase screen (no parameters)
        aps1_obj = Points(file_path=join(self.path, "p1_aps.h5"), logger=self.logger)
        aps1_obj.open(
            other_file_path=join(self.path, "p1_ts_filt.h5"),
            input_path=self.config.general.input_path
        )

        # select second-order points
        cand_mask2 = selectPixels(
            path=self.path, selection_method="temp_coh",
            thrsh=self.config.filtering.coherence_p2,
            grid_size=None, bool_plot=True,
            logger=self.logger
        )  # first-order points are included in second-order points

        slc_stack_obj = slcStack(join(self.config.general.input_path, "slcStack.h5"))
        slc_stack_obj.open()  # for metadata
        pix_sel_obj = BaseStack(file=join(self.path, "selection_method.h5"), logger=self.logger)
        pix_sel_obj.writeToFile(dataset_name="TPC", data=cand_mask2, metadata=slc_stack_obj.metadata, mode="w")
        pix_sel_obj.prepareDataset(dataset_name="TCS", dshape=(point1_obj.length, point1_obj.width), dtype=np.bool_,
                                   metadata=slc_stack_obj.metadata, mode="a")
        pix_sel_obj.prepareDataset(dataset_name="PL", dshape=(point1_obj.length, point1_obj.width), dtype=np.bool_,
                                   metadata=slc_stack_obj.metadata, mode="a")

        if self.config.temporarily_coherent_scatterer.use_tcs:
            self.logger.info(msg="Add Temporarily Coherent Scatterers (TCS).")
            coh_value_tcs = int(self.config.temporarily_coherent_scatterer.coherence_p2_tcs * 100)

            fname = join(self.path, self.config.temporarily_coherent_scatterer.path_tcs,
                         f"lifetime_{self.config.temporarily_coherent_scatterer.method_name}.h5")
            tcs_coh_obj = BaseStack(file=fname, logger=self.logger)
            tcs_coh = tcs_coh_obj.read(dataset_name="coherence_map")
            tcs_length = tcs_coh_obj.read(dataset_name="subset_length_map")

            fig = plt.figure(figsize=(15, 5))
            ax = fig.add_subplot()
            im = ax.imshow(tcs_coh, cmap=cmc.cm.cmaps["grayC"], vmin=0, vmax=1)
            auto_flip_direction(slc_stack_obj.metadata, ax=ax, print_msg=True)
            ax.set_xlabel("Range")
            ax.set_ylabel("Azimuth")
            plt.colorbar(im, pad=0.03, shrink=0.5)
            plt.title("TCS: Temporal coherence")
            plt.tight_layout()
            fig.savefig(join(self.path, "pic", "step_3_tcs_temporal_phase_coherence.png"), dpi=300)
            plt.close(fig)

            fname = join(self.path, self.config.temporarily_coherent_scatterer.path_tcs,
                         f"change_map_{self.config.temporarily_coherent_scatterer.method_name}.h5")
            change_idx_obj = BaseStack(file=fname, logger=self.logger)
            change_idx = change_idx_obj.read(dataset_name="change_index")

            fig = plt.figure(figsize=(15, 5))
            ax = fig.add_subplot()
            cmap = cmc.cm.cmaps["batlow"]
            cmap.set_bad(color='white')  # Set the color for masked values
            masked_data = np.ma.masked_where(change_idx <= 0, change_idx)
            im = ax.imshow(masked_data, vmin=0, vmax=slc_stack_obj.numDate, cmap=cmap, interpolation="nearest")
            auto_flip_direction(slc_stack_obj.metadata, ax=ax, print_msg=True)
            ax.set_xlabel("Range")
            ax.set_ylabel("Azimuth")
            plt.colorbar(im, pad=0.03, shrink=0.5)
            plt.title("TCS: Change time index")
            plt.tight_layout()
            fig.savefig(join(self.path, "pic", "step_3_tcs_change_index.png"), dpi=300)
            plt.close(fig)
            del change_idx, change_idx_obj

            cand_mask_tcs = (tcs_coh > self.config.temporarily_coherent_scatterer.coherence_p2_tcs)
            cand_mask_tcs = cand_mask_tcs & (
                tcs_length >= self.config.temporarily_coherent_scatterer.min_lifetime_length)

            pix_sel_obj.writeToFileBlock(dataset_name="TCS", data=cand_mask_tcs)

            # combine TCS pixels with TPC cand_mask2
            cand_mask2 = cand_mask2 | cand_mask_tcs

            fig = plt.figure(figsize=(15, 5))
            ax = fig.add_subplot()
            bmap_obj.plot(ax=ax, logger=self.logger)
            coord_xy = np.array(np.where(cand_mask_tcs)).transpose()
            val = np.ones_like(cand_mask_tcs)
            sc = ax.scatter(coord_xy[:, 1], coord_xy[:, 0], c=val[cand_mask_tcs], s=0.5,
                            cmap=cmc.cm.cmaps["lajolla_r"],
                            vmin=1, vmax=2)  # set min, max to ensure that points are yellow
            cbar = plt.colorbar(sc, pad=0.03, shrink=0.5)
            cbar.ax.set_visible(False)  # make size of axis consistent with all others
            plt.tight_layout()
            plt.title("Temporarily coherent scatterers")
            fig.savefig(join(self.path, "pic", "step_3_tcs_coh{}.png".format(coh_value_tcs)), dpi=300)
            plt.close(fig)

        if (self.config.phase_linking.use_phase_linking_results &
           (not self.config.temporarily_coherent_scatterer.use_tcs)):
            # read PL results
            pl_coh = readfile.read(join(self.config.phase_linking.inverted_path, "phase_series.h5"),
                                   datasetName='temporalCoherence')[0]
            pl_coh = pl_coh[1, :, :]
            siblings = readfile.read(join(self.config.phase_linking.inverted_path, "phase_series.h5"),
                                     datasetName='shp')[0]

            if self.config.phase_linking.use_ps:
                mask_ps = readfile.read(self.config.phase_linking.mask_ps_file,
                                        datasetName='mask')[0].astype(np.bool_)
                cand_mask_pl = (pl_coh > self.config.filtering.coherence_p2) | mask_ps
            else:
                cand_mask_pl = (pl_coh > self.config.filtering.coherence_p2)
                # remove ps candidates, because the ps detection strategy in miaplpy seems to be biased.
                cand_mask_pl[siblings <= self.config.phase_linking.num_siblings] = False

            if self.config.phase_linking.mask_phase_linking_file is not None:
                path_mask_pl_aoi = join(self.config.phase_linking.mask_phase_linking_file)
                self.logger.info(msg="load mask for area of interest from: {}.".format(path_mask_pl_aoi))
                mask_pl_aoi = readfile.read(path_mask_pl_aoi, datasetName='mask')[0].astype(np.bool_)

                fig = plt.figure(figsize=(15, 5))
                ax = fig.add_subplot()
                ax.imshow(mask_pl_aoi, cmap=cmc.cm.cmaps["grayC"], alpha=0.5, zorder=10, vmin=0, vmax=1)
                bmap_obj.plot(ax=ax, logger=self.logger)
                coord_xy = np.array(np.where(cand_mask_pl)).transpose()
                val = np.ones_like(cand_mask_pl)
                sc = ax.scatter(coord_xy[:, 1], coord_xy[:, 0], c=val[cand_mask_pl], s=0.5,
                                cmap=cmc.cm.cmaps["lajolla_r"],
                                vmin=1, vmax=2)  # set min, max to ensure that points are yellow
                cbar = plt.colorbar(sc, pad=0.03, shrink=0.5)
                cbar.ax.set_visible(False)  # make size of axis consistent with all others
                plt.tight_layout()
                plt.title("Mask for phase linking points")
                fig.savefig(join(self.path, "pic", "step_3_mask_p2_coh{}_phase_linking.png".format(coh_value)), dpi=300)
                plt.close(fig)

                # mask points after plotting, so that available coherent points are visible in figure
                cand_mask_pl[~mask_pl_aoi] = False

            # combine phase linking coherence with TPC cand_mask2
            cand_mask2 = cand_mask2 | cand_mask_pl

            pix_sel_obj.writeToFileBlock(dataset_name="PL", data=cand_mask_pl)

        mask_valid_area = ut.detectValidAreas(bmap_obj=bmap_obj, logger=self.logger)

        pix_sel_obj.writeToFile(dataset_name="valid", data=mask_valid_area, mode="a")

        if self.config.filtering.mask_p2_file is not None:
            path_mask_aoi = join(self.config.filtering.mask_p2_file)
            self.logger.info(msg="load mask for area of interest from: {}.".format(path_mask_aoi))
            mask_aoi = readfile.read(path_mask_aoi, datasetName='mask')[0].astype(np.bool_)
            mask_valid_area &= mask_aoi
            # todo: add unstable points from p1 for densification
        else:
            self.logger.info(msg="No mask for area of interest given.")

        cand_mask2[p1_mask] = True  # add all selected 1.order points to avoid spatial gaps in 2D unwrapping
        # cand_mask2[cand_mask_sparse] = True  # add only stable points from 1.order points

        cand_mask2 &= mask_valid_area

        fig = plt.figure(figsize=(15, 5))
        ax = fig.add_subplot()
        ax.imshow(mask_valid_area, cmap=cmc.cm.cmaps["grayC"], alpha=0.5, zorder=10, vmin=0, vmax=1)
        bmap_obj.plot(ax=ax, logger=self.logger)
        coord_xy = np.array(np.where(cand_mask2)).transpose()
        val = np.ones_like(cand_mask2)
        sc = ax.scatter(coord_xy[:, 1], coord_xy[:, 0], c=val[cand_mask2], s=0.5, cmap=cmc.cm.cmaps["lajolla_r"],
                        vmin=0, vmax=10)  # set min, max to ensure that points are yellow
        cbar = plt.colorbar(sc, pad=0.03, shrink=0.5)
        cbar.ax.set_visible(False)  # make size of axis consistent with all others
        plt.tight_layout()
        plt.title("Mask for dense point set")
        if self.config.temporarily_coherent_scatterer.use_tcs:
            fig_name = f"step_3_mask_p2_coh{coh_value}-{coh_value_tcs}.png"
        else:
            fig_name = f"step_3_mask_p2_coh{coh_value}.png"
        fig.savefig(join(self.path, "pic", fig_name), dpi=300)
        plt.close(fig)

        if self.config.temporarily_coherent_scatterer.use_tcs:
            fname = f"p2_coh{coh_value}-{coh_value_tcs}_ifg_wr.h5"
        else:
            fname = f"p2_coh{coh_value}_ifg_wr.h5"
        point2_obj = Points(file_path=join(self.path, fname), logger=self.logger)
        coord_xy = np.array(np.where(cand_mask2)).transpose()
        point_id2 = point_id_img[cand_mask2]
        point2_obj.prepare(
            point_id=point_id2,
            coord_xy=coord_xy,
            input_path=self.config.general.input_path
        )

        ifg_stack_obj = BaseStack(file=join(self.path, "ifg_stack.h5"), logger=self.logger)

        point2_obj.phase = ut.readPhasePatchwise(stack_obj=ifg_stack_obj, dataset_name="ifgs",
                                                 num_patches=self.config.general.num_patches, cand_mask=cand_mask2,
                                                 point_id_img=point_id_img, logger=self.logger)

        if (self.config.phase_linking.use_phase_linking_results &
           (not self.config.temporarily_coherent_scatterer.use_tcs)):
            self.logger.info(msg="read phase from MiaplPy results...")
            phase_linking_obj = BaseStack(
                file=join(self.config.phase_linking.inverted_path, "phase_series.h5"),
                logger=self.logger
            )

            pl_phase = ut.readPhasePatchwise(
                stack_obj=phase_linking_obj, dataset_name="phase",
                num_patches=self.config.general.num_patches,
                cand_mask=cand_mask2,
                point_id_img=point_id_img, logger=self.logger
            )

            # subset to time span
            slc_stack_obj = slcStack(join(self.config.general.input_path, "slcStack.h5"))
            slc_stack_obj.open()
            time_mask = createTimeMaskFromDates(
                start_date=self.config.preparation.start_date,
                stop_date=self.config.preparation.end_date,
                date_list=slc_stack_obj.dateList,
                logger=self.logger
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

        point2_obj.writeToFile()
        del point2_obj, ifg_stack_obj

        if self.config.temporarily_coherent_scatterer.use_tcs:
            fname = f"p2_coh{coh_value}-{coh_value_tcs}_aps.h5"
            fname_other_file = f"p2_coh{coh_value}-{coh_value_tcs}_ifg_wr.h5"
        else:
            fname = f"p2_coh{coh_value}_aps.h5"
            fname_other_file = f"p2_coh{coh_value}_ifg_wr.h5"
        aps2_obj = Points(file_path=join(self.path, fname), logger=self.logger)
        aps2_obj.open(
            other_file_path=join(self.path, fname_other_file),
            input_path=self.config.general.input_path
        )

        if self.config.filtering.apply_aps_filtering:
            # spatial filtering of points with linear motion only (no non-linear motion)
            if self.config.filtering.interpolation_method == "kriging":
                aps1_phase, aps2_phase = estimateAtmosphericPhaseScreen(
                    residuals=phase_for_aps_filtering,
                    coord_map1=point1_obj.coord_map,
                    coord_map2=aps2_obj.coord_map,
                    num_cores=self.config.general.num_cores,
                    bool_plot=False,
                    logger=self.logger
                )
            else:
                aps1_phase, aps2_phase = simpleInterpolation(
                    residuals=phase_for_aps_filtering,
                    coord_map1=point1_obj.coord_map,
                    coord_map2=aps2_obj.coord_map,
                    interp_method=self.config.filtering.interpolation_method
                )
        else:
            msg = "#" * 10
            msg += " SKIP ATMOSPHERIC FILTERING! "
            msg += "#" * 10
            self.logger.info(msg=msg)
            num_points1 = phase_for_aps_filtering.shape[0]
            num_points2 = aps2_obj.coord_map.shape[0]
            num_time = phase_for_aps_filtering.shape[1]
            aps1_phase = np.zeros((num_points1, num_time), dtype=np.float32)
            aps2_phase = np.zeros((num_points2, num_time), dtype=np.float32)

        point1_obj.phase -= aps1_phase
        point1_obj.writeToFile()

        aps1_obj.phase = aps1_phase
        aps2_obj.phase = aps2_phase
        aps1_obj.writeToFile()
        aps2_obj.writeToFile()

    def runDensificationTimeAndSpace(self):
        """RunDensificationTimeAndSpace."""
        coh_value = int(self.config.filtering.coherence_p2 * 100)

        point2_obj = Points(file_path=join(self.path, "p2_coh{}_ifg_unw.h5".format(coh_value)), logger=self.logger)
        point2_obj.open(
            other_file_path=join(self.path, "p2_coh{}_ifg_wr.h5".format(coh_value)),
            input_path=self.config.general.input_path
        )  # wrapped phase

        # estimate parameters from unwrapped phase
        point1_obj = Points(file_path=join(self.path, "p1_ifg_unw.h5"), logger=self.logger)
        point1_obj.open(input_path=self.config.general.input_path)
        vel_p1, demerr_p1 = ut.estimateParameters(obj=point1_obj, ifg_space=True)[:2]

        # load wrapped phase to remove known components for unwrapping p2 points
        point1_obj = Points(file_path=join(self.path, "p1_ifg_wr.h5"), logger=self.logger)  # wrapped phase!
        point1_obj.open(input_path=self.config.general.input_path)

        aps1_obj = Points(file_path=join(self.path, "p1_aps.h5"), logger=self.logger)
        aps1_obj.open(input_path=self.config.general.input_path)

        aps2_obj = Points(file_path=join(self.path, "p2_coh{}_aps.h5".format(coh_value)), logger=self.logger)
        aps2_obj.open(input_path=self.config.general.input_path)

        if self.config.filtering.mask_p2_file is None:
            """
            overview of points contained in the *_obj
            (unstable p1 means: p1 which were not used in atmospheric filtering)
            p2:   p2 - inconsistent p2 + unstable p1 + stable p1        --> p2:   p2
            aps2: p2 + unstable p1 + stable p1                          --> aps2: p2
            p1:   stable p1 + unstable p1                               --> p1:   stable p1 + unstable p1
            aps1: stable p1                                             --> aps1: stable p1 + unstable p1
            """
            # find unstable p1 in p2 (and in aps2)
            point_id_img = np.arange(0, point1_obj.length * point1_obj.width).reshape(
                (point1_obj.length, point1_obj.width))
            p1_mask = point1_obj.createMask()
            aps1_mask = aps1_obj.createMask()
            mask_unstable_p1 = p1_mask & (~aps1_mask)
            unstable_p1_id = point_id_img[np.where(mask_unstable_p1)]

            mask_unstable_p1_in_p2 = np.ones((aps2_obj.num_points,), dtype=np.bool_)
            for p in aps2_obj.point_id:
                if p not in unstable_p1_id:
                    mask_unstable_p1_in_p2[aps2_obj.point_id == p] = False

            # add unstable p1 from aps2 to aps1
            aps1_obj.addPointsFromObj(
                new_point_id=aps2_obj.point_id[mask_unstable_p1_in_p2],
                new_coord_xy=aps2_obj.coord_xy[mask_unstable_p1_in_p2, :],
                new_phase=aps2_obj.phase[mask_unstable_p1_in_p2, :],
                new_num_points=mask_unstable_p1_in_p2[mask_unstable_p1_in_p2].shape[0],
                input_path=self.config.general.input_path
            )

            # remove unstable p1 from p2 and aps2. thereby remove inconsistent p2 from aps2.
            p2_mask = point2_obj.createMask()
            mask_only_p2 = p2_mask & (~p1_mask)
            keep_id = point_id_img[np.where(mask_only_p2)]
            point2_obj.removePoints(keep_id=keep_id, input_path=self.config.general.input_path)
            aps2_obj.removePoints(keep_id=keep_id, input_path=self.config.general.input_path)

        else:
            """
            if spatial mask is applied:
            p2:   p2 - inconsistent p2 (?) + p1 (coincidently?)         --> p2:   p2
            aps2: p2 + p1 (coincidently?)                               --> aps2: p2
            p1:   stable p1 + unstable p1                               --> p1:   stable p1 (+ unstable p1)
            aps1: stable p1                                             --> aps1: stable p1 (+ unstable p1)
            """
            use_all_p1 = False  # todo: add to config
            if use_all_p1:
                raise NotImplementedError("Use all p1 is not implemented.")
            else:
                # remove also values from estimated parameters
                mask = np.ones((point1_obj.num_points,), dtype=np.bool_)
                for p in point1_obj.point_id:
                    if p not in aps1_obj.point_id:
                        mask[point1_obj.point_id == p] = False

                vel_p1 = vel_p1[mask]
                demerr_p1 = demerr_p1[mask]

                # remove unstable p1 from p1
                point1_obj.removePoints(
                    keep_id=aps1_obj.point_id,
                    input_path=self.config.general.input_path
                )

                # remove p2 which are coincidentally equal to p1
                point_id_img = np.arange(0, point1_obj.length * point1_obj.width).reshape(
                    (point1_obj.length, point1_obj.width))
                p1_mask = point1_obj.createMask()
                p2_mask = point2_obj.createMask()
                mask_p2 = ~(p1_mask & p2_mask) & p2_mask
                p2_id = point_id_img[np.where(mask_p2)]
                point2_obj.removePoints(keep_id=p2_id, input_path=self.config.general.input_path)
                aps2_obj.removePoints(keep_id=p2_id, input_path=self.config.general.input_path)

        # return to ifg-space
        a_ifg = point2_obj.ifg_net_obj.getDesignMatrix()
        aps1_ifg_phase = np.matmul(a_ifg, aps1_obj.phase.T).T
        aps2_ifg_phase = np.matmul(a_ifg, aps2_obj.phase.T).T

        # correct for APS
        point2_obj.phase = np.angle(np.exp(1j * point2_obj.phase) * np.conjugate(np.exp(1j * aps2_ifg_phase)))
        point1_obj.phase = np.angle(np.exp(1j * point1_obj.phase) * np.conjugate(np.exp(1j * aps1_ifg_phase)))

        demerr, vel, gamma = densifyNetwork(
            point1_obj=point1_obj,
            vel_p1=vel_p1,
            demerr_p1=demerr_p1,
            point2_obj=point2_obj,
            num_conn_p1=self.config.densification.num_connections_to_p1,
            max_dist_p1=self.config.densification.max_distance_to_p1,
            velocity_bound=self.config.densification.velocity_bound,
            demerr_bound=self.config.densification.dem_error_bound,
            num_samples=self.config.densification.num_optimization_samples,
            num_cores=self.config.general.num_cores,
            logger=self.logger
        )  # returns parameters of both first- and second-order points

        # store combined set of first and second-order points
        point2_obj.addPointsFromObj(
            new_point_id=point1_obj.point_id,
            new_coord_xy=point1_obj.coord_xy,
            new_phase=point1_obj.phase,
            new_num_points=point1_obj.num_points,
            input_path=self.config.general.input_path
        )

        bmap_obj = AmplitudeImage(file_path=join(self.path, "background_map.h5"))
        fig = viewer.plotScatter(value=gamma, coord=point2_obj.coord_xy, bmap_obj=bmap_obj,
                                 ttl="Coherence from temporal unwrapping\nBefore outlier removal", s=3.5,
                                 cmap="lajolla", vmin=0, vmax=1, logger=self.logger)[0]
        fig.savefig(join(self.path, "pic", "step_4_temporal_unwrapping_p2_coh{}.png".format(coh_value)), dpi=300)
        plt.close(fig)

        mask_gamma = gamma >= self.config.densification.arc_unwrapping_coherence
        self.logger.info(msg=f"Reduce the dense point set by {mask_gamma[~mask_gamma].shape[0]} points,")
        self.logger.info(msg=f"due to coherence from temporal unwrapping < "
                             f"{self.config.densification.arc_unwrapping_coherence}")
        point2_obj.removePoints(mask=mask_gamma, keep_id=[], input_path=self.config.general.input_path)

        fig = plt.figure(figsize=(15, 5))
        axs = fig.subplots(1, 2)
        axs[0].hist(-vel[mask_gamma] * 100, bins=200)
        axs[0].set_ylabel('Absolute frequency')
        axs[0].set_xlabel('Mean velocity [cm / year]')

        axs[1].hist(-demerr[mask_gamma], bins=200)
        axs[1].set_ylabel('Absolute frequency')
        axs[1].set_xlabel('DEM error [m]')
        fig.savefig(join(self.path, "pic", "step_4_consistency_parameters_p2_coh{}.png".format(coh_value)),
                    dpi=300)
        plt.close(fig)

        fig = viewer.plotScatter(value=gamma[mask_gamma], coord=point2_obj.coord_xy, bmap_obj=bmap_obj,
                                 ttl="Coherence from temporal unwrapping\nAfter outlier removal", s=3.5,
                                 cmap="lajolla", vmin=0, vmax=1, logger=self.logger)[0]
        fig.savefig(join(self.path, "pic", "step_4_temporal_unwrapping_p2_coh{}_reduced.png".format(coh_value)),
                    dpi=300)
        plt.close(fig)

        fig = viewer.plotScatter(value=-vel[mask_gamma], coord=point2_obj.coord_xy,
                                 ttl="Mean velocity in [m / year]",
                                 bmap_obj=bmap_obj, s=3.5, cmap="roma", symmetric=True,
                                 logger=self.logger)[0]
        fig.savefig(join(self.path, "pic", "step_4_estimation_velocity_p2_coh{}.png".format(coh_value)), dpi=300)
        plt.close(fig)

        fig = viewer.plotScatter(value=-demerr[mask_gamma], coord=point2_obj.coord_xy, ttl="DEM correction in [m]",
                                 bmap_obj=bmap_obj, s=3.5, cmap="vanimo", symmetric=True,
                                 logger=self.logger)[0]
        fig.savefig(join(self.path, "pic", "step_4_estimation_dem_correction_p2_coh{}.png".format(coh_value)), dpi=300)
        plt.close(fig)

        self.logger.info(msg="Remove phase contributions from mean velocity "
                             "and DEM correction from wrapped phase of points.")
        pred_phase_demerr, pred_phase_vel = ut.predictPhase(
            obj=point2_obj,
            vel=vel[mask_gamma],
            demerr=demerr[mask_gamma],
            ifg_space=True,
            logger=self.logger
        )
        pred_phase = pred_phase_demerr + pred_phase_vel

        wr_phase = point2_obj.phase
        wr_res_phase = np.angle(np.exp(1j * wr_phase) * np.conjugate(np.exp(1j * pred_phase)))

        triang_obj = PointNetworkTriangulation(coord_xy=point2_obj.coord_xy, coord_map_xy=None, logger=self.logger)
        triang_obj.triangulateGlobal()
        arcs = triang_obj.getArcsFromAdjMat()

        unw_res_phase = spatialUnwrapping(num_ifgs=point2_obj.ifg_net_obj.num_ifgs,
                                          num_points=point2_obj.num_points,
                                          phase=wr_res_phase,
                                          method=self.config.general.spatial_unwrapping_method,
                                          edges=arcs,
                                          num_cores=self.config.general.num_cores, logger=self.logger)

        self.logger.info(msg="Add phase contributions from mean velocity "
                             "and DEM correction back to spatially unwrapped residual phase.")
        unw_phase = unw_res_phase + pred_phase

        point2_obj.phase = unw_phase
        vel = ut.estimateParameters(obj=point2_obj, ifg_space=True)[0]
        point2_obj.phase = ut.setReferenceToPeakOfHistogram(phase=unw_phase, vel=vel, num_bins=300)

        point2_obj.writeToFile()

        phase_ts = ut.invertIfgNetwork(
            phase=unw_phase,
            num_points=point2_obj.num_points,
            ifg_net_obj=point2_obj.ifg_net_obj,
            logger=self.logger)

        point_obj = Points(file_path=join(self.path, "p2_coh{}_ts.h5".format(coh_value)), logger=self.logger)
        point_obj.open(
            other_file_path=join(self.path, "p2_coh{}_ifg_unw.h5".format(coh_value)),
            input_path=self.config.general.input_path
        )
        point_obj.phase = phase_ts

        point_obj.writeToFile()

    def runDensificationSpace(self):
        """RunDensification."""
        coh_value = int(self.config.filtering.coherence_p2 * 100)

        point_obj = Points(file_path=join(self.path, "p2_coh{}_ifg_unw.h5".format(coh_value)), logger=self.logger)
        point_obj.open(
            other_file_path=join(self.path, "p2_coh{}_ifg_wr.h5".format(coh_value)),
            input_path=self.config.general.input_path
        )  # open wr phase

        aps2_obj = Points(file_path=join(self.path, "p2_coh{}_aps.h5".format(coh_value)), logger=self.logger)
        aps2_obj.open(input_path=self.config.general.input_path)

        # return to ifg-space
        a_ifg = point_obj.ifg_net_obj.getDesignMatrix()
        aps2_ifg_phase = np.matmul(a_ifg, aps2_obj.phase.T).T

        # correct for APS
        point_obj.phase = np.angle(np.exp(1j * point_obj.phase) * np.conjugate(np.exp(1j * aps2_ifg_phase)))

        triang_obj = PointNetworkTriangulation(coord_xy=point_obj.coord_xy, coord_map_xy=None, logger=self.logger)
        triang_obj.triangulateGlobal()  # if coord_map is not given, only global delaunay and knn can be calculated
        arcs = triang_obj.getArcsFromAdjMat()

        unw_phase = spatialUnwrapping(num_ifgs=point_obj.ifg_net_obj.num_ifgs,
                                      num_points=point_obj.num_points,
                                      phase=point_obj.phase,
                                      method=self.config.general.spatial_unwrapping_method,
                                      edges=arcs,
                                      num_cores=self.config.general.num_cores, logger=self.logger)

        # adjust reference to peak of histogram
        point_obj.phase = unw_phase
        vel = ut.estimateParameters(obj=point_obj, ifg_space=True)[0]
        point_obj.phase = ut.setReferenceToPeakOfHistogram(phase=unw_phase, vel=vel, num_bins=300)

        point_obj.writeToFile()
        del point_obj

        point_obj = Points(file_path=join(self.path, "p2_coh{}_ts.h5".format(coh_value)), logger=self.logger)
        point_obj.open(
            other_file_path=join(self.path, "p2_coh{}_ifg_wr.h5".format(coh_value)),
            input_path=self.config.general.input_path
        )

        phase_ts = ut.invertIfgNetwork(phase=unw_phase, num_points=point_obj.num_points,
                                       ifg_net_obj=point_obj.ifg_net_obj,
                                       logger=self.logger)

        point_obj.phase = phase_ts

        point_obj.writeToFile()
