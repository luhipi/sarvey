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

"""Tests for SARvey `ifg_network` functions."""

import os
import unittest
from datetime import datetime
import sys
import logging

import numpy as np

from sarvey.ifg_network import SmallBaselineNetwork


class TestUtils(unittest.TestCase):
    root_path = None
    config_file = None
    configuration = None
    output_data_path = None

    @classmethod
    def setUp(cls) -> None:
        """Define the Class method SetUp."""
        # define class variables, create temporary output dir etc. here
        cls.root_path = "./"
        if os.path.basename(os.getcwd()) == "tests":
            cls.root_path = "../"

        cls.logger = logging.getLogger(__name__)
        console_handler = logging.StreamHandler(sys.stdout)
        cls.logger.addHandler(console_handler)
        cls.logger.setLevel(logging.getLevelName('DEBUG'))

    @classmethod
    def tearDown(cls) -> None:
        """Define the Class method tearDown."""
        # delete testfolder or do some other cleanup here

    def testConfigure_ok(self):
        """Test the expected output."""
        # Input:
        pbase = np.array([0, 0, 0, 0])  # not important for this test
        dates = [datetime(2023, 8, 17), datetime(2023, 8, 17), datetime(2023, 8, 17), datetime(2023, 8, 17)]

        tbase = np.array([0, 6, 12, 18])
        ifg_net_obj = SmallBaselineNetwork(logger=self.logger)
        ifg_net_obj.configure(pbase=pbase, tbase=tbase, num_link=2, max_tbase=12, dates=dates)
        ifg_list = np.array([(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)])
        assert (ifg_net_obj.ifg_list == ifg_list).all()

        tbase = np.array([0, 12, 18, 312])
        ifg_net_obj = SmallBaselineNetwork(logger=self.logger)
        ifg_net_obj.configure(pbase=pbase, tbase=tbase, num_link=3, max_tbase=5, dates=dates)
        ifg_list = np.array([(0, 1), (1, 2), (2, 3)])
        assert (ifg_net_obj.ifg_list == ifg_list).all()

        tbase = np.array([0, 12, 18, 312])
        ifg_net_obj = SmallBaselineNetwork(logger=self.logger)
        ifg_net_obj.configure(pbase=pbase, tbase=tbase, num_link=3, max_tbase=20, dates=dates)
        ifg_list = np.array([(0, 1), (0, 2), (1, 2), (2, 3)])
        assert (ifg_net_obj.ifg_list == ifg_list).all()

    # def testConfigure_err(self):
    #     """Test for expected Errors."""
    #
