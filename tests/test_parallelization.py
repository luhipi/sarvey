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

"""Tests for worker-count utilities used by SARvey parallel processing."""

import os
import unittest
import numpy as np

from sarvey.utils import clampNumCores, scaleNumCores, splitDatasetForParallelProcessing, setNativeThreadEnvIfUnset


class TestParallelizationUtils(unittest.TestCase):
    """Test safety guards for parallel worker allocation."""

    _native_thread_vars = (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    )

    def setUp(self):
        self._env_backup = {k: os.environ.get(k) for k in self._native_thread_vars}

    def tearDown(self):
        for key, val in self._env_backup.items():
            if val is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = val

    def testClampNumCores_caps_to_num_samples(self):
        """Clamp worker count to available samples."""
        self.assertEqual(clampNumCores(requested_cores=16, num_samples=3), 3)
        self.assertEqual(clampNumCores(requested_cores=4, num_samples=10), 4)

    def testClampNumCores_rejects_invalid_values(self):
        """Reject non-positive workers or samples."""
        with self.assertRaises(ValueError):
            clampNumCores(requested_cores=0, num_samples=5)
        with self.assertRaises(ValueError):
            clampNumCores(requested_cores=2, num_samples=0)

    def testScaleNumCores_keeps_at_least_one_worker(self):
        """Avoid the old floor-to-zero behavior for requested cores below ten."""
        self.assertEqual(scaleNumCores(requested_cores=2, num_samples=50, scale=0.1), 1)
        self.assertEqual(scaleNumCores(requested_cores=9, num_samples=50, scale=0.1), 1)
        self.assertEqual(scaleNumCores(requested_cores=10, num_samples=50, scale=0.1), 1)
        self.assertEqual(scaleNumCores(requested_cores=20, num_samples=50, scale=0.1), 2)

    def testSplitDatasetForParallelProcessing_clamps_num_cores(self):
        """Split logic should automatically cap workers to sample count."""
        idx = splitDatasetForParallelProcessing(num_samples=3, num_cores=8)
        self.assertEqual(len(idx), 3)
        np.testing.assert_array_equal(idx[0], np.array([0]))
        np.testing.assert_array_equal(idx[1], np.array([1]))
        np.testing.assert_array_equal(idx[2], np.array([2]))

    def testSplitDatasetForParallelProcessing_balances_ranges(self):
        """Chunks should remain balanced when sample count is not divisible by workers."""
        idx = splitDatasetForParallelProcessing(num_samples=10, num_cores=3)
        lengths = [cur.shape[0] for cur in idx]
        self.assertEqual(lengths, [4, 3, 3])
        np.testing.assert_array_equal(np.concatenate(idx), np.arange(10))

    def testSetNativeThreadEnvIfUnset_sets_only_missing(self):
        """Unset variables are initialized while pre-set variables are kept."""
        os.environ.pop("OMP_NUM_THREADS", None)
        os.environ["OPENBLAS_NUM_THREADS"] = "7"
        os.environ.pop("MKL_NUM_THREADS", None)
        os.environ.pop("NUMEXPR_NUM_THREADS", None)

        updates = setNativeThreadEnvIfUnset(num_threads=1)

        self.assertEqual(os.environ["OMP_NUM_THREADS"], "1")
        self.assertEqual(os.environ["OPENBLAS_NUM_THREADS"], "7")
        self.assertEqual(os.environ["MKL_NUM_THREADS"], "1")
        self.assertEqual(os.environ["NUMEXPR_NUM_THREADS"], "1")
        self.assertEqual(updates["OMP_NUM_THREADS"], "1")
        self.assertEqual(updates["MKL_NUM_THREADS"], "1")
        self.assertEqual(updates["NUMEXPR_NUM_THREADS"], "1")
        self.assertTrue("OPENBLAS_NUM_THREADS" not in updates)

    def testSetNativeThreadEnvIfUnset_rejects_invalid_threads(self):
        """Native thread limit must be positive."""
        with self.assertRaises(ValueError):
            setNativeThreadEnvIfUnset(num_threads=0)


if __name__ == '__main__':
    unittest.main()
