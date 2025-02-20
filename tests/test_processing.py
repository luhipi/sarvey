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

"""Tests for `SARvey` package."""
import logging
import os
import shutil
import sys
import unittest
from copy import deepcopy
from glob import glob
from os.path import join

from sarvey.sarvey_mti import createParser, run
from sarvey.config import loadConfiguration


class TestProcessing(unittest.TestCase):
    root_path = None
    config_file = None
    configuration = None
    output_data_path = None
    args = None
    logger = None

    @classmethod
    def setUp(self) -> None:
        """Define the Class method SetUp."""
        # define class variables, create temporary output dir etc. here
        self.root_path = "./"
        if os.path.basename(os.getcwd()) == "tests":
            self.root_path = "../"

        self.config_file = os.path.abspath(f"{self.root_path}/tests/testdata/config_test.json")
        self.configuration = loadConfiguration(path=self.config_file)

        self.logger = logging.getLogger(__name__)
        console_handler = logging.StreamHandler(sys.stdout)
        self.logger.addHandler(console_handler)
        self.logger.setLevel(logging.getLevelName('DEBUG'))

        self.output_data_path = os.path.abspath(f"{self.root_path}/tests/testdata/output")
        parser = createParser()
        self.args = parser.parse_args(["-f", self.config_file, "0", "4", "-w", self.output_data_path])

        if not os.path.exists(self.output_data_path):
            os.makedirs(join(self.output_data_path, "pic"))
        else:
            if not os.path.exists(join(self.output_data_path, "pic")):
                os.mkdir(join(self.output_data_path, "pic"))

    @classmethod
    def tearDown(self) -> None:
        """Define the Class method tearDown."""
        # delete testfolder or do some other cleanup here
        try:
            if os.path.exists(self.output_data_path):
                shutil.rmtree(self.output_data_path)
        except OSError:
            print("Deletion of the directory %s failed" % self.output_data_path)
        else:
            print("Successfully deleted the directory %s" % self.output_data_path)

    def testUnwrappingTimeAndSpace(self):
        """TestUnwrappingTimeAndSpace."""
        config = deepcopy(self.configuration)
        args = deepcopy(self.args)

        # set config
        config.general.apply_temporal_unwrapping = True
        config.filtering.apply_aps_filtering = False

        # preparation
        args.start = 0
        args.stop = 0

        run(config=config, args=args, logger=self.logger)
        assert glob(os.path.join(config.general.output_path, "background_map.h5")), \
            'Processing failed (background_map.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "coordinates_map.h5")), \
            'Processing failed (coordinates_map.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "ifg_network.h5")), \
            'Processing failed (ifg_network.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "ifg_stack.h5")), \
            'Processing failed (ifg_stack.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "temporal_coherence.h5")), \
            'Processing failed (temporal_coherence.h5 not created).'

        # consistencyCheck
        args.start = 1
        args.stop = 1

        run(config=config, args=args, logger=self.logger)
        assert glob(os.path.join(config.general.output_path, "point_network.h5")), \
            'Processing failed (point_network.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "point_network_parameter.h5")), \
            'Processing failed (point_network_parameter.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "p1_ifg_wr.h5")), \
            'Processing failed (p1_ifg_wr.h5 not created).'

        # unwrapping
        args.start = 2
        args.stop = 2

        run(config=config, args=args, logger=self.logger)
        assert glob(os.path.join(config.general.output_path, "p1_ifg_unw.h5")), \
            'Processing failed (p1_ifg_unw.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "p1_ts.h5")), \
            'Processing failed (p1_ts.h5 not created).'

        # filtering
        args.start = 3
        args.stop = 3

        run(config=config, args=args, logger=self.logger)
        assert glob(os.path.join(config.general.output_path, "p1_ts_filt.h5")), \
            'Processing failed (p1_ts_filt.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "p1_aps.h5")), \
            'Processing failed (p1_aps.h5 not created).'
        coh_value = int(config.filtering.coherence_p2 * 100)
        assert glob(os.path.join(config.general.output_path, f"p2_coh{coh_value}_ifg_wr.h5")), \
            f'Processing failed (p2_coh{coh_value}_ifg_wr.h5 not created).'
        assert glob(os.path.join(config.general.output_path, f"p2_coh{coh_value}_aps.h5")), \
            f'Processing failed (p2_coh{coh_value}_aps.h5 not created).'

        # densification
        args.start = 4
        args.stop = 4

        run(config=config, args=args, logger=self.logger)
        coh_value = int(config.filtering.coherence_p2 * 100)
        assert glob(os.path.join(config.general.output_path, f"p2_coh{coh_value}_ifg_unw.h5")), \
            f'Processing failed (p2_coh{coh_value}_ifg_unw.h5 not created).'
        assert glob(os.path.join(config.general.output_path, f"p2_coh{coh_value}_ts.h5")), \
            f'Processing failed (p2_coh{coh_value}_ts.h5 not created).'

    def testUnwrappingSpace(self):
        """TestUnwrappingTimeAndSpace."""
        config = deepcopy(self.configuration)
        args = deepcopy(self.args)

        # set config
        config.general.apply_temporal_unwrapping = False
        config.preparation.ifg_network_type = "sb"

        # preparation
        args.start = 0
        args.stop = 0

        run(config=config, args=args, logger=self.logger)
        assert glob(os.path.join(config.general.output_path, "background_map.h5")), \
            'Processing failed (background_map.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "coordinates_map.h5")), \
            'Processing failed (coordinates_map.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "ifg_network.h5")), \
            'Processing failed (ifg_network.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "ifg_stack.h5")), \
            'Processing failed (ifg_stack.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "temporal_coherence.h5")), \
            'Processing failed (temporal_coherence.h5 not created).'

        # consistencyCheck
        args.start = 1
        args.stop = 1

        run(config=config, args=args, logger=self.logger)
        assert glob(os.path.join(config.general.output_path, "point_network.h5")), \
            'Processing failed (point_network.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "point_network_parameter.h5")), \
            'Processing failed (point_network_parameter.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "p1_ifg_wr.h5")), \
            'Processing failed (p1_ifg_wr.h5 not created).'

        # unwrapping
        args.start = 2
        args.stop = 2

        run(config=config, args=args, logger=self.logger)
        assert glob(os.path.join(config.general.output_path, "p1_ifg_unw.h5")), \
            'Processing failed (p1_ifg_unw.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "p1_ts.h5")), \
            'Processing failed (p1_ts.h5 not created).'

        # filtering
        args.start = 3
        args.stop = 3

        run(config=config, args=args, logger=self.logger)
        assert glob(os.path.join(config.general.output_path, "p1_ts_filt.h5")), \
            'Processing failed (p1_ts_filt.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "p1_aps.h5")), \
            'Processing failed (p1_aps.h5 not created).'
        coh_value = int(config.filtering.coherence_p2 * 100)
        assert glob(os.path.join(config.general.output_path, f"p2_coh{coh_value}_ifg_wr.h5")), \
            f'Processing failed (p2_coh{coh_value}_ifg_wr.h5 not created).'
        assert glob(os.path.join(config.general.output_path, f"p2_coh{coh_value}_aps.h5")), \
            f'Processing failed (p2_coh{coh_value}_aps.h5 not created).'

        # densification
        args.start = 4
        args.stop = 4

        run(config=config, args=args, logger=self.logger)
        coh_value = int(config.filtering.coherence_p2 * 100)
        assert glob(os.path.join(config.general.output_path, f"p2_coh{coh_value}_ifg_unw.h5")), \
            f'Processing failed (p2_coh{coh_value}_ifg_unw.h5 not created).'
        assert glob(os.path.join(config.general.output_path, f"p2_coh{coh_value}_ts.h5")), \
            f'Processing failed (p2_coh{coh_value}_ts.h5 not created).'

    def testPhaseLinking(self):
        """TestUnwrappingTimeAndSpace."""
        config = deepcopy(self.configuration)
        args = deepcopy(self.args)

        # set config
        config.phase_linking.use_phase_linking_results = True
        config.phase_linking.use_ps = True
        config.general.apply_temporal_unwrapping = True
        config.filtering.coherence_p2 = 0.75

        # preparation
        args.start = 0
        args.stop = 0

        run(config=config, args=args, logger=self.logger)
        assert glob(os.path.join(config.general.output_path, "background_map.h5")), \
            'Processing failed (background_map.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "coordinates_map.h5")), \
            'Processing failed (coordinates_map.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "ifg_network.h5")), \
            'Processing failed (ifg_network.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "ifg_stack.h5")), \
            'Processing failed (ifg_stack.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "temporal_coherence.h5")), \
            'Processing failed (temporal_coherence.h5 not created).'

        # consistencyCheck
        args.start = 1
        args.stop = 1

        run(config=config, args=args, logger=self.logger)
        assert glob(os.path.join(config.general.output_path, "point_network.h5")), \
            'Processing failed (point_network.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "point_network_parameter.h5")), \
            'Processing failed (point_network_parameter.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "p1_ifg_wr.h5")), \
            'Processing failed (p1_ifg_wr.h5 not created).'

        # unwrapping
        args.start = 2
        args.stop = 2

        run(config=config, args=args, logger=self.logger)
        assert glob(os.path.join(config.general.output_path, "p1_ifg_unw.h5")), \
            'Processing failed (p1_ifg_unw.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "p1_ts.h5")), \
            'Processing failed (p1_ts.h5 not created).'

        # filtering
        args.start = 3
        args.stop = 3

        run(config=config, args=args, logger=self.logger)
        assert glob(os.path.join(config.general.output_path, "p1_ts_filt.h5")), \
            'Processing failed (p1_ts_filt.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "p1_aps.h5")), \
            'Processing failed (p1_aps.h5 not created).'
        coh_value = int(config.filtering.coherence_p2 * 100)
        assert glob(os.path.join(config.general.output_path, f"p2_coh{coh_value}_ifg_wr.h5")), \
            f'Processing failed (p2_coh{coh_value}_ifg_wr.h5 not created).'
        assert glob(os.path.join(config.general.output_path, f"p2_coh{coh_value}_aps.h5")), \
            f'Processing failed (p2_coh{coh_value}_aps.h5 not created).'

        # densification
        args.start = 4
        args.stop = 4

        run(config=config, args=args, logger=self.logger)
        coh_value = int(config.filtering.coherence_p2 * 100)
        assert glob(os.path.join(config.general.output_path, f"p2_coh{coh_value}_ifg_unw.h5")), \
            f'Processing failed (p2_coh{coh_value}_ifg_unw.h5 not created).'
        assert glob(os.path.join(config.general.output_path, f"p2_coh{coh_value}_ts.h5")), \
            f'Processing failed (p2_coh{coh_value}_ts.h5 not created).'

    def testMasking(self):
        """TestUnwrappingTimeAndSpace."""
        config = deepcopy(self.configuration)
        args = deepcopy(self.args)

        # set config
        config.general.apply_temporal_unwrapping = False
        config.preparation.ifg_network_type = "sb"
        config.consistency_check.mask_p1_file = "tests/testdata/aoi_mask.h5"
        config.filtering.mask_p2_file = "tests/testdata/aoi_mask.h5"

        # preparation
        args.start = 0
        args.stop = 0

        run(config=config, args=args, logger=self.logger)
        assert glob(os.path.join(config.general.output_path, "background_map.h5")), \
            'Processing failed (background_map.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "coordinates_map.h5")), \
            'Processing failed (coordinates_map.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "ifg_network.h5")), \
            'Processing failed (ifg_network.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "ifg_stack.h5")), \
            'Processing failed (ifg_stack.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "temporal_coherence.h5")), \
            'Processing failed (temporal_coherence.h5 not created).'

        # consistencyCheck
        args.start = 1
        args.stop = 1

        run(config=config, args=args, logger=self.logger)
        assert glob(os.path.join(config.general.output_path, "point_network.h5")), \
            'Processing failed (point_network.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "point_network_parameter.h5")), \
            'Processing failed (point_network_parameter.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "p1_ifg_wr.h5")), \
            'Processing failed (p1_ifg_wr.h5 not created).'

        # unwrapping
        args.start = 2
        args.stop = 2

        run(config=config, args=args, logger=self.logger)
        assert glob(os.path.join(config.general.output_path, "p1_ifg_unw.h5")), \
            'Processing failed (p1_ifg_unw.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "p1_ts.h5")), \
            'Processing failed (p1_ts.h5 not created).'

        # filtering
        args.start = 3
        args.stop = 3

        run(config=config, args=args, logger=self.logger)
        assert glob(os.path.join(config.general.output_path, "p1_ts_filt.h5")), \
            'Processing failed (p1_ts_filt.h5 not created).'
        assert glob(os.path.join(config.general.output_path, "p1_aps.h5")), \
            'Processing failed (p1_aps.h5 not created).'
        coh_value = int(config.filtering.coherence_p2 * 100)
        assert glob(os.path.join(config.general.output_path, f"p2_coh{coh_value}_ifg_wr.h5")), \
            f'Processing failed (p2_coh{coh_value}_ifg_wr.h5 not created).'
        assert glob(os.path.join(config.general.output_path, f"p2_coh{coh_value}_aps.h5")), \
            f'Processing failed (p2_coh{coh_value}_aps.h5 not created).'

        # densification
        args.start = 4
        args.stop = 4

        run(config=config, args=args, logger=self.logger)
        coh_value = int(config.filtering.coherence_p2 * 100)
        assert glob(os.path.join(config.general.output_path, f"p2_coh{coh_value}_ifg_unw.h5")), \
            f'Processing failed (p2_coh{coh_value}_ifg_unw.h5 not created).'
        assert glob(os.path.join(config.general.output_path, f"p2_coh{coh_value}_ts.h5")), \
            f'Processing failed (p2_coh{coh_value}_ts.h5 not created).'
