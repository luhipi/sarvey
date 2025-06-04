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

"""IfgNetwork module for SARvey."""
import datetime
import h5py
import os
import matplotlib.pyplot as plt
import numpy as np
from typing import Union
import warnings
from logging import Logger
from scipy.spatial import Delaunay


class IfgNetwork:
    """Abstract class/interface for different types of interferogram networks."""

    ifg_list: Union[list, np.ndarray] = None

    def __init__(self):
        """Init."""
        self.pbase = None
        self.tbase = None
        self.num_images = None
        self.ifg_list = list()  # is later converted to np.array
        self.pbase_ifg = None
        self.tbase_ifg = None
        self.num_ifgs = None
        self.dates = list()

    def plot(self):
        """Plot the network of interferograms."""
        fig = plt.figure(figsize=(15, 5))
        axs = fig.subplots(1, 3)
        dt = [datetime.date.fromisoformat(d) for d in self.dates]
        axs[0].plot(dt, self.pbase, 'ko')
        for idx in self.ifg_list:
            xx = np.array([dt[idx[0]], dt[idx[1]]])
            yy = np.array([self.pbase[idx[0]], self.pbase[idx[1]]])
            axs[0].plot(xx, yy, 'k-')
        axs[0].set_ylabel('perpendicular baseline [m]')
        axs[0].set_xlabel('temporal baseline [years]')
        axs[0].set_title('Network of interferograms')
        fig.autofmt_xdate()

        axs[1].hist(self.tbase_ifg * 365.25, bins=100)
        axs[1].set_ylabel('Absolute frequency')
        axs[1].set_xlabel('temporal baseline [days]')

        axs[2].hist(self.pbase_ifg, bins=100)
        axs[2].set_ylabel('Absolute frequency')
        axs[2].set_xlabel('perpendicular baseline [m]')
        return fig

    def getDesignMatrix(self):
        """Compute the design matrix for the smallbaseline network."""
        a = np.zeros((self.num_ifgs, self.num_images))
        for i in range(len(self.ifg_list)):
            a[i, self.ifg_list[i][0]] = 1
            a[i, self.ifg_list[i][1]] = -1
        return a

    def open(self, *, path: str):
        """Read stored information from already existing.h5 file.

        Parameter
        -----------
        path: str
            path to existing file to read from.
        """
        with h5py.File(path, 'r') as f:
            self.num_images = f.attrs["num_images"]
            self.num_ifgs = f.attrs["num_ifgs"]

            self.tbase_ifg = f['tbase_ifg'][:]
            self.pbase_ifg = f['pbase_ifg'][:]
            self.tbase = f['tbase'][:]
            self.pbase = f['pbase'][:]
            self.ifg_list = f['ifg_list'][:]
            try:
                self.dates = f['dates'][:]
                self.dates = [date.decode("utf-8") for date in self.dates]
            except KeyError as ke:
                self.dates = None
                print(f"IfgNetwork is in old dataformat. Cannot read 'dates'! {ke}")

            f.close()

    def writeToFile(self, *, path: str, logger: Logger):
        """Write all existing data to .h5 file.

        Parameters
        ----------
        path: str
            path to filename
        logger: Logger
            Logging handler.
        """
        logger.info(msg="write IfgNetwork to {}".format(path))

        if os.path.exists(path):
            os.remove(path)

        dates = np.array(self.dates, dtype=np.bytes_)

        with h5py.File(path, 'w') as f:
            f.attrs["num_images"] = self.num_images
            f.attrs["num_ifgs"] = self.num_ifgs

            f.create_dataset('tbase_ifg', data=self.tbase_ifg)
            f.create_dataset('pbase_ifg', data=self.pbase_ifg)
            f.create_dataset('tbase', data=self.tbase)
            f.create_dataset('pbase', data=self.pbase)
            f.create_dataset('ifg_list', data=self.ifg_list)
            f.create_dataset('dates', data=dates)


class StarNetwork(IfgNetwork):
    """Star network of interferograms (single-reference)."""

    def configure(self, *, pbase: np.ndarray, tbase: np.ndarray, ref_idx: int, dates: list):
        """Create list of interferograms containing the indices of the images and computes baselines.

        Parameter
        ---------
        pbase: np.ndarray
            Perpendicular baselines of the SAR acquisitions.
        tbase: np.ndarray
            Temporal baselines of the SAR acquisitions.
        ref_idx: int
            Index of the reference image.
        dates: list
            Dates of the acquisitions.
        """
        self.pbase = pbase
        self.tbase = tbase / 365.25
        self.num_images = pbase.shape[0]
        self.dates = dates

        for i in range(self.num_images):
            if i == ref_idx:
                continue
            self.ifg_list.append((ref_idx, i))

        self.pbase_ifg = np.delete(self.pbase - self.pbase[ref_idx], ref_idx)
        self.tbase_ifg = np.delete(self.tbase - self.tbase[ref_idx], ref_idx)
        self.num_ifgs = self.num_images - 1


class SmallTemporalBaselinesNetwork(IfgNetwork):
    """Small temporal baselines network of interferograms without restrictions on the perpendicular baselines."""

    def configure(self, *, pbase: np.ndarray, tbase: np.ndarray, num_link: int = None, dates: list):
        """Create list of interferograms containing the indices of the images and computes baselines.

        Parameter
        -----------
        pbase: np.ndarray
            Perpendicular baselines of the SAR acquisitions.
        tbase: np.ndarray
            Temporal baselines of the SAR acquisitions.
        num_link: int
            Number of consecutive links in time connecting acquisitions.
        dates: list
            Dates of the acquisitions.
        """
        self.pbase = pbase
        self.tbase = tbase / 365.25
        self.num_images = pbase.shape[0]
        self.dates = dates

        for i in range(self.num_images):
            for j in range(num_link):
                if i + j + 1 >= self.num_images:
                    continue
                self.ifg_list.append((i, i + j + 1))

        self.ifg_list = [(i, j) for i, j in self.ifg_list if i != j]  # remove connections to itself, e.g. (0, 0)

        self.pbase_ifg = np.array([self.pbase[idx[1]] - self.pbase[idx[0]] for idx in self.ifg_list])
        self.tbase_ifg = np.array([self.tbase[idx[1]] - self.tbase[idx[0]] for idx in self.ifg_list])
        self.num_ifgs = self.pbase_ifg.shape[0]


class SmallBaselineNetwork(IfgNetwork):
    """Small baseline network of interferograms restricting both temporal and spatial baselines."""

    def configure(self, *, pbase: np.ndarray, tbase: np.ndarray, num_link: int, max_tbase: int, dates: list):
        """Create list of interferograms containing the indices of the images and computes baselines.

        Parameter
        -----------
        pbase: np.ndarray
            perpendicular baselines of the SAR acquisitions.
        tbase: np.ndarray
            temporal baselines of the SAR acquisitions.
        max_tbase: int
            maximum temporal baseline in [days] (default: None).
        num_link: int
            number of links within the range of maximum temporal baseline.
        dates: list
            Dates of the acquisitions.
        """
        self.pbase = pbase
        self.tbase = tbase / 365.25
        self.num_images = pbase.shape[0]
        self.dates = dates
        flag_restrict_to_max_tbase = False

        # in this section use tbase in [days] (function argument, not self.)
        for i in range(self.num_images - 1):

            if i + 1 < self.num_images - 1:
                # always use one connection to nearest neighbour in time
                self.ifg_list.append((i, i + 1))
            else:
                self.ifg_list.append((i, i + 1))
                break
            # compute index corresponding to max_tbase for current time
            diff = np.abs(tbase - (tbase[i] + max_tbase))
            max_idx = np.where(diff == diff.min())[0][0]
            self.ifg_list.append((i, max_idx))

            if max_idx == i:  # no further images between i and max_idx
                flag_restrict_to_max_tbase = True
                continue

            # spread the rest of the links over the remaining time steps in between
            links = np.floor(np.arange(i, max_idx, (max_idx - i) / (num_link - 1)))[1:].astype(int)
            for link in links:
                self.ifg_list.append((i, link))
        self.ifg_list = np.unique(self.ifg_list, axis=0)

        if flag_restrict_to_max_tbase:
            warnings.warn(f"Cannot restrict ifgs to maximum temporal baseline of {max_tbase} days.")

        self.ifg_list = [(i, j) for i, j in self.ifg_list if i != j]  # remove connections to itself, e.g. (0, 0)

        self.pbase_ifg = np.array([self.pbase[idx[1]] - self.pbase[idx[0]] for idx in self.ifg_list])
        self.tbase_ifg = np.array([self.tbase[idx[1]] - self.tbase[idx[0]] for idx in self.ifg_list])
        self.num_ifgs = self.pbase_ifg.shape[0]


class DelaunayNetwork(IfgNetwork):
    """Delaunay network of interferograms which restricts both the temporal and perpendicular baselines."""

    def configure(self, *, pbase: np.ndarray, tbase: np.ndarray, dates: list):
        """Create list of interferograms containing the indices of the images and computes baselines.

        Parameter
        -----------
        pbase: np.ndarray
            perpendicular baselines of the SAR acquisitions, array
        tbase: np.ndarray
            temporal baselines of the SAR acquisitions, array
        dates: list
            Dates of the acquisitions, list.
        """
        self.pbase = pbase
        self.tbase = tbase / 365.25
        self.num_images = pbase.shape[0]
        self.dates = dates
        scale = 0.25

        network = Delaunay(points=np.stack([self.pbase, self.tbase * 365.25 * scale]).T)
        for p1, p2, p3 in network.simplices:
            self.ifg_list.append((p1, p2))
            self.ifg_list.append((p1, p3))
            self.ifg_list.append((p2, p3))
        self.ifg_list = np.unique(self.ifg_list, axis=0)

        self.pbase_ifg = np.array([self.pbase[idx[1]] - self.pbase[idx[0]] for idx in self.ifg_list])
        self.tbase_ifg = np.array([self.tbase[idx[1]] - self.tbase[idx[0]] for idx in self.ifg_list])
        self.num_ifgs = self.pbase_ifg.shape[0]


class SmallBaselineYearlyNetwork(IfgNetwork):
    """Small baseline network of interferograms with yearly connections."""

    def configure(self, *, pbase: np.ndarray, tbase: np.ndarray, num_link: int = None, dates: list):
        """Create list of interferograms containing the indices of the images and computes baselines.

        Parameter
        -----------
        pbase: np.ndarray
            perpendicular baselines of the SAR acquisitions, array
        tbase: np.ndarray
            temporal baselines of the SAR acquisitions, array
        num_link: int
            Number of consecutive links in time connecting acquisitions.
        dates: list
            Dates of the acquisitions, list.
        """
        self.pbase = pbase
        self.tbase = tbase / 365.25
        self.num_images = pbase.shape[0]
        self.dates = dates

        # add small temporal baselines
        for i in range(self.num_images):
            for j in range(num_link):
                if i + j + 1 >= self.num_images:
                    continue
                self.ifg_list.append((i, i + j + 1))

        # add yearly ifgs
        for i in range(self.num_images):
            # find index of image at roughly one year distance
            diff = np.abs(tbase - (tbase[i] + 365.25))
            year_idx = np.where(diff == diff.min())[0][0]
            print(year_idx)
            if year_idx != self.num_images - 1:  # avoid connections to the last image
                self.ifg_list.append((i, year_idx))
                print("found!")

        self.ifg_list = np.unique(self.ifg_list, axis=0)
        self.ifg_list = [(i, j) for i, j in self.ifg_list if i != j]  # remove connections to itself, e.g. (0, 0)

        self.pbase_ifg = np.array([self.pbase[idx[1]] - self.pbase[idx[0]] for idx in self.ifg_list])
        self.tbase_ifg = np.array([self.tbase[idx[1]] - self.tbase[idx[0]] for idx in self.ifg_list])
        self.num_ifgs = self.pbase_ifg.shape[0]
