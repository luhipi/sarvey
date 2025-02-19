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

"""Preparation module for SARvey."""
import datetime
import matplotlib.pyplot as plt
import numpy as np
from logging import Logger
from os.path import join

import mintpy.utils.readfile as readfile

from sarvey import viewer
import sarvey.utils as ut
from sarvey.objects import CoordinatesUTM, AmplitudeImage, BaseStack, Points
from sarvey.triangulation import PointNetworkTriangulation


def createTimeMaskFromDates(*, start_date: str, stop_date: str, date_list: list, logger: Logger):
    """Create a mask with selected dates within given time frame.

    Parameters
    ----------
    start_date: str
        Start date.
    stop_date: str
        Stop date.
    date_list: list
        all avaiable dates in the slcStack.h5.
    logger: Logger
        Logging handler.

    Returns
    -------
    time_mask: np.ndarray
        mask with True for selected dates.
    num_slc: int
        number of selected images.
    result_date_list: list
        list of selected dates.
    """
    time_mask = np.ones((len(date_list)), dtype=np.bool_)
    date_list = [datetime.date(year=int(d[:4]), month=int(d[4:6]), day=int(d[6:])) for d in date_list]

    if (start_date is None) and (stop_date is None):
        # use all images.
        result_date_list = [date.isoformat() for date in date_list]
        return time_mask, time_mask.shape[0], result_date_list

    if start_date is None:
        start_date = min(date_list)
    else:
        start_date = datetime.date.fromisoformat(start_date)

    if stop_date is None:
        stop_date = max(date_list)
    else:
        stop_date = datetime.date.fromisoformat(stop_date)

    if start_date >= stop_date:
        logger.error(msg="Choose start date < stop date!")
        raise ValueError

    if stop_date < min(date_list):
        logger.error(msg="Stop date is before the first acquired image. Choose a later stop date!")
        raise ValueError

    if start_date > max(date_list):
        logger.error(msg="Start date is after the last acquired image. Choose an earlier start date!")
        raise ValueError

    shift = "    "
    logger.debug(msg=shift + "{:>10} {:>10}".format("   Date   ", "Selected"))
    logger.debug(msg=shift + "{:>10} {:>10}".format("__________", "________"))

    result_date_list = list()
    for i, date in enumerate(date_list):
        if (date < start_date) or (date > stop_date):
            time_mask[i] = False
        else:
            result_date_list.append(date.isoformat())
        val = "      x" if time_mask[i] else ""
        logger.debug(msg=shift + "{:>10} {:>3}".format(date.isoformat(), val))

    num_slc = time_mask[time_mask].shape[0]
    return time_mask, num_slc, result_date_list


def readSlcFromMiaplpy(*, path: str, box: tuple = None, logger: Logger) -> np.ndarray:
    """Read SLC data from phase-linking results of Miaplpy.

    Parameters
    ----------
    path: str
        Path to the phase_series.h5 file.
    box: tuple
        Bounding Box to read from.
    logger: Logger
        Logging handler.

    Returns
    -------
    slc: np.ndarray
        slc stack created from phase-linking results.
    """
    logger.info(msg="read phase from MiaplPy results...")
    phase = readfile.read(path, datasetName='phase', box=box)[0]

    logger.info(msg="read amplitude from MiaplPy results...")
    amp = readfile.read(path, datasetName='amplitude', box=box)[0]

    logger.info(msg="combine phase and amplitude to slc...")
    slc = amp * np.exp(phase * 1j)
    return slc


def readCoherenceFromMiaplpy(*, path: str, box: tuple = None, logger: Logger) -> tuple[np.ndarray, dict]:
    """Read the coherence image from phase-linking of MiaplPy.

    Parameters
    ----------
    path: str
        Path to phase_series.h5 file.
    box: tuple
        Bounding Box to read from.
    logger: Logger
        Logging handler.

    Returns
    -------
    temp_coh: np.ndarray
        temporal coherence image from phase-linking results of MiaplPy.
    """
    logger.info(msg="read quality from MiaplPy results...")
    temp_coh = readfile.read(path, datasetName='temporalCoherence', box=box)[0][1, :, :]
    return temp_coh


def selectPixels(*, path: str, selection_method: str, thrsh: float,
                 grid_size: int = None, bool_plot: bool = False, logger: Logger):
    """Select pixels based on temporal coherence.

    Parameters
    ----------
    path: str
        Path to the directory with the temporal_coherence.h5 file.
    selection_method: str
        Pixel selection method. Currently, only "temp_coh" is implemented.
    thrsh: float
        Threshold for pixel selection.
    grid_size: int
        Grid size for sparse pixel selection.
    bool_plot: bool
        Plot the selected pixels.
    logger: Logger
        Logging handler.

    Returns
    -------
    cand_mask: np.ndarray
        Mask with selected pixels.
    """
    quality = None
    grid_min_val = None
    cand_mask = None
    unit = None
    cmap = None
    # compute candidates
    if selection_method == "temp_coh":
        temp_coh_obj = BaseStack(file=join(path, "temporal_coherence.h5"), logger=logger)
        quality = temp_coh_obj.read(dataset_name="temp_coh")
        cand_mask = quality >= thrsh
        grid_min_val = False
        unit = "Temporal\nCoherence [ ]"
        cmap = "lajolla"

    if selection_method == "miaplpy":
        raise NotImplementedError("This part is not developed yet. MiaplPy data is read in another way.")
        # pl_coherence = readCoherenceFromMiaplpy(path=join(path, 'inverted', 'phase_series.h5'), box=None,
        # logger=logger)
        # cand_mask = pl_coherence >= thrsh
        # quality = pl_coherence
        # grid_min_val = False
        # unit = "Phase-Linking\nCoherence [ ]"
        # cmap = "lajolla"

    if grid_size is not None:  # -> sparse pixel selection
        coord_utm_obj = CoordinatesUTM(file_path=join(path, "coordinates_utm.h5"), logger=logger)
        coord_utm_obj.open()
        box_list = ut.createSpatialGrid(coord_utm_img=coord_utm_obj.coord_utm,
                                        length=coord_utm_obj.coord_utm.shape[1],
                                        width=coord_utm_obj.coord_utm.shape[2],
                                        grid_size=grid_size,
                                        logger=logger)[0]
        cand_mask_sparse = ut.selectBestPointsInGrid(box_list=box_list, quality=quality, sel_min=grid_min_val)
        cand_mask &= cand_mask_sparse

    if bool_plot:
        coord_xy = np.array(np.where(cand_mask)).transpose()
        bmap_obj = AmplitudeImage(file_path=join(path, "background_map.h5"))
        viewer.plotScatter(value=quality[cand_mask], coord=coord_xy, bmap_obj=bmap_obj, ttl="Selected pixels",
                           unit=unit, s=2, cmap=cmap, vmin=0, vmax=1, logger=logger)
        # if grid_size is not None:
        #     psViewer.plotGridFromBoxList(box_list, ax=ax, edgecolor="k", linewidth=0.2)
        plt.tight_layout()
        plt.gcf().savefig(join(path, "pic", "selected_pixels_{}_{}.png".format(selection_method, thrsh)),
                          dpi=300)
        plt.close(plt.gcf())

    return cand_mask


def createArcsBetweenPoints(*, point_obj: Points, knn: int = None, max_arc_length: float = np.inf,
                            logger: Logger) -> np.ndarray:
    """Create a spatial network of arcs to triangulate the points.

    All points are triangulated with a Delaunay triangulation. If knn is given, the triangulation is done with the k
    nearest neighbors. Too long arcs are removed from the network. If, afterward, the network is not connected, a
    delaunay triangulation is performed again to ensure connectivity in the network.

    Parameters
    ----------
    point_obj: Points
        Point object.
    knn: int
        Number of nearest neighbors to consider (default: None).
    max_arc_length: float
        Maximum length of an arc. Longer arcs will be removed. Default: np.inf.
    logger: Logger
        Logging handler.

    Returns
    -------
    arcs: np.ndarray
        Arcs of the triangulation containing the indices of the points for each arc.
    """
    triang_obj = PointNetworkTriangulation(coord_xy=point_obj.coord_xy, coord_utmxy=point_obj.coord_utm, logger=logger)

    if knn is not None:
        triang_obj.triangulateKnn(k=knn)

    triang_obj.triangulateGlobal()

    logger.info(msg="remove arcs with length > {}.".format(max_arc_length))
    triang_obj.removeLongArcs(max_dist=max_arc_length)

    if not triang_obj.isConnected():
        triang_obj.triangulateGlobal()

    logger.info(msg="retrieve arcs from adjacency matrix.")
    arcs = triang_obj.getArcsFromAdjMat()
    return arcs
