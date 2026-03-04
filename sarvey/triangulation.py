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

"""Triangulation module for SARvey."""
import time
import numpy as np
from scipy.spatial import Delaunay, KDTree
from scipy.sparse import lil_matrix
from logging import Logger

from mintpy.utils import ptime


class PointNetworkTriangulation:
    """PointNetworkTriangulation."""

    def __init__(self, *, coord_xy: np.ndarray, logger: Logger):
        """Triangulate points in space based on distance.

        Parameters
        ----------
        coord_xy: np.ndarray
            Radar coordinates of the points.
        logger: Logger
            Logging handler.
        """
        self.coord_xy = coord_xy
        num_points = self.coord_xy.shape[0]
        self.logger = logger

        self.logger.debug(f"Initializing Point Network Triangulation with {num_points} points...")

        # create sparse matrix with dim (num_points x num_points), add 1 if connected.
        # at the end create the network once to reduce computational time.
        self.adj_mat = lil_matrix((num_points, num_points), dtype=np.bool_)

    def getArcsFromAdjMat(self):
        """Convert the adjacency matrix into a list of arcs.

        Returns
        -------
        arcs: np.ndarray
            List of arcs with indices of the start and end point.
        """
        self.logger.debug(f"Extracting arcs from adjacency matrix with size {self.adj_mat.shape}...")
        a = self.adj_mat.copy()
        # copy entries from lower to upper triangular matrix
        b = (a + a.T)
        # remove entries from diagonal and lower part of matrix
        arc_tmp = [[i, b.indices[b.indptr[i]:b.indptr[i + 1]]] for i in range(b.shape[0])]
        arc_tmp = [[s, e_list[np.where(e_list < s)[0]]] for s, e_list in arc_tmp]

        arcs = list()
        for s, e_list in arc_tmp:
            for e in e_list:
                arcs.append([s, e])
        arcs = np.array(arcs)
        return arcs

    def triangulateDelaunay(self):
        """Connect the points with a delaunay triangulation."""
        self.logger.debug("Triangulating points with delaunay...")

        network = Delaunay(points=self.coord_xy)
        self.logger.debug(f"Number of triangles in Delaunay triangulation: {network.simplices.shape[0]}")
        for p1, p2, p3 in network.simplices:
            self.adj_mat[p1, p2] = True
            self.adj_mat[p1, p3] = True
            self.adj_mat[p2, p3] = True

    def triangulateCircularNearestNeighbors(self, *, num_partitions: int = 8):
        """Connect points to the nearest neighbors in all directions based on circular partitions.

        Parameters
        ----------
        num_partitions: int
            Number of partitions to divide the area around each point. Default: 8.
        """
        self.logger.info(msg="Triangulate points with circular nearest neighbors.")
        num_points = self.coord_xy.shape[0]
        angles = np.linspace(0, 2 * np.pi, num_partitions + 1)[:-1]
        tree = KDTree(data=self.coord_xy)

        prog_bar = ptime.progressBar(maxValue=num_points)
        start_time = time.time()

        k_search = min(150, num_points)
        for p1 in range(num_points):
            # find all nearest neighbours independent of the direction within the max_dist
            idx = tree.query(self.coord_xy[p1, :], k=k_search)[1]
            idx = np.atleast_1d(idx)
            idx = idx[1:]

            # add the neighbours to the predefined bins based on the angle to the point
            angles_to_neighbors = np.arctan2(self.coord_xy[idx, 1] - self.coord_xy[p1, 1],
                                             self.coord_xy[idx, 0] - self.coord_xy[p1, 0])
            angle_bins = np.digitize(angles_to_neighbors, angles) - 1  # -1 to make it zero-indexed

            # select the closest neighbour per angle partition
            for angle_idx in range(num_partitions):
                # find the indices of the neighbours in this angle partition
                angle_mask = (angle_bins == angle_idx)
                if np.any(angle_mask):
                    # take the neighbour from the bin which has the smallest idx number
                    closest_idx = np.min(idx[angle_mask])

                    # add the connection to the adjacency matrix
                    self.adj_mat[p1, closest_idx] = True
            prog_bar.update(value=p1 + 1, every=np.int16(num_points / (num_points / 5)),
                            suffix='{}/{} points triangulated'.format(p1 + 1, num_points + 1))
        prog_bar.close()
        m, s = divmod(time.time() - start_time, 60)
        self.logger.info(f"time used for triangulation: {m:02.0f} mins {s:02.1f} secs.")
