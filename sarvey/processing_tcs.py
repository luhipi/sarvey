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
"""Temporarily Coherent Scatterer processing module for SARvey."""
from os.path import join
import numpy as np
from logging import Logger

from sarvey.objects import Points, BaseStack
from sarvey.tcs_utils import selectCoherentImagesAndIfgs, preTriangulateEachIfg
from sarvey.unwrapping import spatialUnwrapping
import sarvey.utils as ut
from sarvey.config import Config


def runDensificationSpace(*, path: str, config: Config, logger: Logger):
    """RunDensification."""
    coh_value = int(config.filtering.coherence_p2 * 100)
    coh_value_tcs = int(config.temporarily_coherent_scatterer.coherence_p2_tcs * 100)

    point_obj = Points(file_path=join(path, f"p2_coh{coh_value}-{coh_value_tcs}_ifg_unw.h5"), logger=logger)
    point_obj.open(
        other_file_path=join(path, f"p2_coh{coh_value}-{coh_value_tcs}_ifg_wr.h5"),
        input_path=config.general.input_path
    )  # open wr phase

    aps2_obj = Points(file_path=join(path, f"p2_coh{coh_value}-{coh_value_tcs}_aps.h5"), logger=logger)
    aps2_obj.open(input_path=config.general.input_path)

    # return to ifg-space
    a_ifg = point_obj.ifg_net_obj.getDesignMatrix()
    aps2_ifg_phase = np.matmul(a_ifg, aps2_obj.phase.T).T

    # correct for APS
    point_obj.phase = np.angle(np.exp(1j * point_obj.phase) * np.conjugate(np.exp(1j * aps2_ifg_phase)))

    # identify coherent lifetime of TCS
    fname = join(config.temporarily_coherent_scatterer.tcs_path,
                 f"lifetime_{config.temporarily_coherent_scatterer.method_name}.h5")
    tcs_coh_obj = BaseStack(file=fname, logger=logger)
    subset_index_map = tcs_coh_obj.read(dataset_name="subset_index_map")

    fname = join(config.temporarily_coherent_scatterer.tcs_path,
                 f"change_map_{config.temporarily_coherent_scatterer.method_name}.h5")
    tcs_change_idx_obj = BaseStack(file=fname,
                                   logger=logger)
    change_index_map = tcs_change_idx_obj.read(dataset_name="change_index")

    sca_type_obj = BaseStack(file=join(path, "scatterer_type.h5"), logger=logger)
    mask_ccs = sca_type_obj.read(dataset_name="CCS")
    mask_tcs = sca_type_obj.read(dataset_name="TCS")

    lifetime_images, lifetime_ifgs = selectCoherentImagesAndIfgs(
        change_index_map=change_index_map,
        subset_index_map=subset_index_map,
        point_obj=point_obj,
        mask_ccs=mask_ccs,
        mask_tcs=mask_tcs
    )

    edges_per_ifg = preTriangulateEachIfg(
        coord_xy=point_obj.coord_xy,
        lifetime_ifgs=lifetime_ifgs,
        num_cores=config.general.num_cores,
        logger=logger
    )

    unw_phase = spatialUnwrapping(
        num_ifgs=point_obj.ifg_net_obj.num_ifgs,
        num_points=point_obj.num_points,
        phase=point_obj.phase,
        lifetime_ifgs=lifetime_ifgs,
        edges=edges_per_ifg,
        method=config.general.spatial_unwrapping_method,
        num_cores=config.general.num_cores,
        logger=logger
    )

    # adjust reference to peak of histogram
    point_obj.phase = unw_phase
    vel = ut.estimateParameters(obj=point_obj, ifg_space=True)[0]
    point_obj.phase = ut.setReferenceToPeakOfHistogram(phase=unw_phase, vel=vel, num_bins=300)

    point_obj.writeToFile()
    del point_obj

    point_obj = Points(file_path=join(path, f"p2_coh{coh_value}-{coh_value_tcs}_ts.h5"), logger=logger)
    point_obj.open(
        other_file_path=join(path, f"p2_coh{coh_value}-{coh_value_tcs}_ifg_unw.h5"),
        input_path=config.general.input_path
    )

    phase_ts = ut.invertIfgNetwork(
        phase=unw_phase,
        lifetime_ifgs=lifetime_ifgs,
        lifetime_images=lifetime_images,
        num_points=point_obj.num_points,
        ifg_net_obj=point_obj.ifg_net_obj,
        logger=logger
    )

    point_obj.phase = phase_ts

    point_obj.writeToFile()
