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

"""Triangulation module for SARvey."""
import time
from typing import Optional
import numpy as np
from scipy.spatial import Delaunay, distance_matrix, KDTree
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.csgraph import connected_components
from logging import Logger

from mintpy.utils import ptime


class PointNetworkTriangulation:
    """PointNetworkTriangulation."""

    def __init__(self, *, coord_xy: np.ndarray, coord_utmxy: Optional[np.ndarray], logger: Logger):
        """Triangulate points in space based on distance.

        Parameters
        ----------
        coord_xy: np.ndarray
            Radar coordinates of the points.
        coord_utmxy: np.ndarray
            UTM coordinates of the points.
        logger: Logger
            Logging handler.
        """
        self.coord_xy = coord_xy
        num_points = self.coord_xy.shape[0]
        self.logger = logger

        # create sparse matrix with dim (num_points x num_points), add 1 if connected.
        # create network afterwards once. reduces time.
        self.adj_mat = lil_matrix((num_points, num_points), dtype=np.bool_)

        if coord_utmxy is not None:
            logger.info("create distance matrix between all points...")
            self.dist_mat = distance_matrix(coord_utmxy, coord_utmxy)
            # todo: check out alternatives:
            #       scipy.spatial.KDTree.sparse_distance_matrix
        else:  # if only global delaunay shall be computed without memory issues
            self.dist_mat = None

    def getArcsFromAdjMat(self):
        """Convert the adjacency matrix into a list of arcs.

        Returns
        -------
        arcs: np.ndarray
            List of arcs with indices of the start and end point.
        """
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

    def removeLongArcs(self, *, max_dist: float):
        """Remove arcs from network which are longer than given threshold.

        Parameter
        ---------
        max_dist: float
            distance threshold on arc length in [m]
        """
        mask = self.dist_mat > max_dist
        self.adj_mat[mask] = False

    def isConnected(self):
        """Check if the network is connected."""
        n_components = connected_components(csgraph=csr_matrix(self.adj_mat), directed=False, return_labels=False)
        if n_components == 1:
            return True
        else:
            return False

    def triangulateGlobal(self):
        """Connect the points with a GLOBAL delaunay triangulation."""
        self.logger.info("Triangulate points with global delaunay.")

        network = Delaunay(points=self.coord_xy)
        for p1, p2, p3 in network.simplices:
            self.adj_mat[p1, p2] = True
            self.adj_mat[p1, p3] = True
            self.adj_mat[p2, p3] = True

    def triangulateKnn(self, *, k: int):
        """Connect points to the k-nearest neighbours."""
        self.logger.info(msg="Triangulate points with {}-nearest neighbours.".format(k))
        num_points = self.coord_xy.shape[0]
        prog_bar = ptime.progressBar(maxValue=num_points)
        start_time = time.time()
        count = 0
        tree = KDTree(data=self.coord_xy)

        if k > num_points:
            k = num_points
            self.logger.info("k > number of points. Connect all points with each other.")
        for p1 in range(num_points):
            idx = tree.query(self.coord_xy[p1, :], k)[1]
            self.adj_mat[p1, idx] = True
            count += 1
            prog_bar.update(value=count + 1, every=np.int16(num_points / 250),
                            suffix='{}/{} points triangulated'.format(count + 1, num_points + 1))
        prog_bar.close()
        m, s = divmod(time.time() - start_time, 60)
        self.logger.debug('time used: {:02.0f} mins {:02.1f} secs.'.format(m, s))
