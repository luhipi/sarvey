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

"""Osm utils module for SARvey."""
import numpy as np
import overpy
from logging import Logger
from shapely import Point

from mintpy.utils import readfile, utils as ut


def getSpatialExtend(*, geom_file: str, logger: Logger):
    """Get spatial extend of the radar image.

    Parameters
    ----------
    geom_file: str
        path of geometryRadar.h5 file
    logger: Logger
        Logging handler.

    Returns
    -------
    ll_corner_wgs: list
        list of coordinates of the lower-left corner of the radar image in WGS84 coordinates.
    ur_corner_wgs: list
        list of coordinates of the upper-right corner of the radar image in WGS84 coordinates.
    coord: np.ndarray
        coordinates of all pixels in the radar image in WGS84.
    atr: dict
        metadata dictionary from geometryRadar.h5.
    """
    logger.info('read spatial extend from geometryRadar.h5')
    _, atr = readfile.read(geom_file)
    coord = ut.coordinate(atr, lookup_file=geom_file)
    lat, atr = readfile.read(geom_file, datasetName='latitude')
    lon, _ = readfile.read(geom_file, datasetName='longitude')

    # radar image is flipped up-down
    # unclear: check if bounding box fits. Otherwise, change to max and min values of lat and lon
    ll_bbox = [np.nanmin(lat), np.nanmin(lon)]
    ur_bbox = [np.nanmax(lat), np.nanmax(lon)]

    img_ext = [
        Point(lon[0, 0], lat[0, 0]),
        Point(lon[-1, 0], lat[-1, 0]),
        Point(lon[-1, -1], lat[-1, -1]),
        Point(lon[0, -1], lat[0, -1])
    ]
    return ll_bbox, ur_bbox, img_ext, coord, atr


def runOsmQuery(*, ll_corner_wgs: np.ndarray, ur_corner_wgs: np.ndarray, type_list: list,
                logger: Logger) -> overpy.Result:
    """Query OSM database for transport infrastructure within the spatial extent of the radar image.

    Parameters
    ----------
    ll_corner_wgs: np.ndarray
        coordinates of the lower-left corner of the radar image in WGS84 coordinates.
    ur_corner_wgs: np.ndarray
        coordinates of the upper-right corner of the radar image in WGS84 coordinates.
    type_list: list
        List of street types that shall be queried at the OSM database.
    logger: Logger
        Logging handler.

    Returns
    -------
        result: overpy.Result
            results of the overpy query to OSM database.
    """
    # Initialize overpass connection
    api = overpy.Overpass()

    # Request data from API
    logger.info('querying OSM database for infra types...')
    # query_cmd = "way({},{},{},{}) [""highway=motorway_link""]; (._;>;); out body;"

    query_cmd = "[bbox: {},{},{},{}];("
    for infra_type in type_list:
        logger.info(f"\t - {infra_type}")
        if infra_type == 'rail':
            query_cmd += "way[railway={}];".format(infra_type)
        else:
            query_cmd += "way[highway={}];".format(infra_type)

    query_cmd += ");(._; >;); out body;"  # skel

    cmd = query_cmd.format(ll_corner_wgs[0], ll_corner_wgs[1],
                           ur_corner_wgs[0], ur_corner_wgs[1])
    logger.info("\n" + cmd + "\n")
    result = api.query(cmd)

    if len(result.ways) == 0:
        logger.error('Empty OSM query results. No roads or railway tracks found.')
        raise ValueError

    logger.info('...done.')
    return result


def runOsmQueryBridge(*, ll_corner_wgs: np.ndarray, ur_corner_wgs: np.ndarray, bridge_highway: bool,
                      bridge_railway: bool, logger: Logger) -> overpy.Result:
    """Query OSM database for bridges of transport infrastructure within the spatial extent of the radar image.

    Parameters
    ----------
    ll_corner_wgs: np.ndarray
        coordinates of the lower-left corner of the radar image in WGS84 coordinates.
    ur_corner_wgs: np.ndarray
        coordinates of the upper-right corner of the radar image in WGS84 coordinates.
    bridge_highway: bool
        Set true to query highway bridges.
    bridge_railway: bool
        Set true to query railway bridges.
    logger: Logger
        Logging handler.

    Returns
    -------
        result: overpy.Result
            results of the overpy query to OSM database.
    """
    # Initialize overpass connection
    api = overpy.Overpass()

    # Request data from API
    logger.info('querying OSM database for infra types...')
    # query_cmd = "way({},{},{},{}) [""highway=motorway_link""]; (._;>;); out body;"

    query_cmd = "[bbox: {},{},{},{}];("

    if bridge_highway:
        logger.info('\t - bridge_highway')
        query_cmd += 'way[highway~"^(motorway|motorway_link|trunk|trunk_link)$"][bridge];'

    if bridge_railway:
        logger.info('\t - bridge_railway')
        query_cmd += 'way[railway=rail][bridge];'

    if (bridge_highway is False) & (bridge_railway is False):
        logger.info('\t - all bridges')
        query_cmd += 'way[bridge];'

    query_cmd += ");(._; >;); out body;"  # skel

    cmd = query_cmd.format(ll_corner_wgs[0], ll_corner_wgs[1],
                           ur_corner_wgs[0], ur_corner_wgs[1])
    logger.info("\n" + cmd + "\n")
    result = api.query(cmd)

    if len(result.ways) == 0:
        logger.error('Empty OSM query results. No bridges found.')
        raise ValueError

    logger.info('...done.')
    return result
