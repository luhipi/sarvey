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

"""Download openstreetmap data for area of interest."""
import argparse
import logging
import os
import sys
import time
from os.path import join
import geopandas as gpd
from shapely import ops, Point
import matplotlib

from sarvey import version
from sarvey.osm_utils import runOsmQueryBridge, runOsmQuery, getSpatialExtend

try:
    matplotlib.use('TkAgg')
except ImportError as e:
    print(e)


EXAMPLE = """Example:
  sarvey_osm --geom ./geometryRadar.h5 --railway                       # download railway
  sarvey_osm --geom ./geometryRadar.h5 --highway                       # download highway
  sarvey_osm --geom ./geometryRadar.h5 --railway --bridge              # download railway bridge
  sarvey_osm --geom ./geometryRadar.h5 --railway -o mask_railway.shp   # specify output path
"""


def create_parser():
    """Create_parser."""
    parser = argparse.ArgumentParser(
        description='Download transport infrastructure information from openstreetmap and store as shp-file.',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=EXAMPLE)

    parser.add_argument('-w', '--work_dir', dest='work_dir', default=None,
                        help='absolute path to working directory\n' +
                             '(default: current directory).')

    parser.add_argument('--geom', dest='geom_file', default=None,
                        help='path to existing geometryRadar.h5 file')

    parser.add_argument('--railway', dest='railway', action="store_true", default=False,
                        help='Set true to query railways.')

    parser.add_argument('--highway', dest='highway', action="store_true", default=False,
                        help='Set true to query highways.')

    parser.add_argument('--bridge', dest='bridge', action="store_true", default=False,
                        help='Set true to mask bridges.\n' +
                             'If --railway or --highway set true, only railway/highway bridges are queried.')

    parser.add_argument('-o', dest='out_file_name', default='osm_infra.shp',
                        help="name of output file. (default: 'osm_infra.shp')")

    parser.add_argument('--version', action='version',
                        version=f"SARvey version {version.__version__} - {version.__versionalias__}, "
                                f"{version.__versiondate__}")

    return parser


def downloadOSM(*, railway: bool, highway: bool, bridge: bool,
                work_dir: str, out_file_name: str, logger: logging.Logger, geom_file: str):
    """Download openstreetmap data and store to file.

    Parameters
    ----------
    railway: bool
        download railway data.
    highway: bool
        download highway data.
    bridge: bool
        download bridge data.
    work_dir: str
        working directory.
    out_file_name: str
        output file name.
    logger: logging.Logger
        logger.
    geom_file: str
        path to geometryRadar.h5 file.
    """
    logger.info(msg="Start creating mask file based on openstreetmap data.")

    # get bounding box
    ll_bbox, ur_bbox, img_ext, coord, atr = getSpatialExtend(geom_file=geom_file, logger=logger)

    # store image extend
    gdf = gpd.GeoDataFrame({"geometry": gpd.geoseries.GeoSeries(img_ext)})
    gdf = gdf.dissolve().convex_hull
    gdf.to_file(join(work_dir, "img_extend.gpkg"))

    # store bounding box
    bbox_points = [
        Point(ll_bbox[1], ll_bbox[0]),
        Point(ur_bbox[1], ll_bbox[0]),
        Point(ur_bbox[1], ur_bbox[0]),
        Point(ll_bbox[1], ur_bbox[0])
    ]

    gdf = gpd.GeoDataFrame({"geometry": gpd.geoseries.GeoSeries(bbox_points)})
    gdf = gdf.dissolve().convex_hull
    gdf.to_file(join(work_dir, "bounding_box.gpkg"))

    if (not railway) & (not highway) & (not bridge):
        logger.error(msg="No infrastructure type was specified.")
        return

    if bridge:
        # get requested OSM layer
        query_result = runOsmQueryBridge(
            ll_corner_wgs=ll_bbox, ur_corner_wgs=ur_bbox,
            bridge_highway=highway, bridge_railway=railway,
            logger=logger
        )
    else:
        type_list = list()
        if railway:
            type_list += ["rail"]
        if highway:
            type_list += ["motorway", "motorway_link", "trunk", "trunk_link"]

        # get requested OSM layer
        query_result = runOsmQuery(ll_corner_wgs=ll_bbox, ur_corner_wgs=ur_bbox,
                                   type_list=type_list, logger=logger)

    multi_line_list = list()
    for way in query_result.ways:
        if "area" in way.tags:
            if way.tags["area"] == "yes":
                logger.info('Area flag is true')
                continue
        else:
            # keep coordinates in lat/lon. It will be needed in masking step.
            coord = [[float(way.nodes[i].lon), float(way.nodes[i].lat)] for i in range(len(way.nodes))]
            multi_line_list.append(coord)

    # Merge all road segments
    merged_road = list(ops.linemerge(multi_line_list).geoms)
    gdf = gpd.GeoDataFrame({"geometry": gpd.GeoSeries(merged_road)})
    # gdf = gdf.set_crs(crs=utm_crs)  # set appropriate CRS
    # todo: add attributes if required

    # todo: check ending of output file name
    gdf.to_file(join(work_dir, out_file_name))
    logger.info(msg="OSM download finished.")


def main(iargs=None):
    """Download openstreetmap data and store to file."""
    # check input
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # initiate logger
    logging_level = logging.getLevelName('DEBUG')

    log_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)

    current_datetime = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
    log_filename = f"sarvey_osm_{current_datetime}.log"
    if not os.path.exists(os.path.join(os.getcwd(), "logfiles")):
        os.mkdir(os.path.join(os.getcwd(), "logfiles"))
    file_handler = logging.FileHandler(filename=os.path.join(os.getcwd(), "logfiles", log_filename))
    file_handler.setFormatter(log_format)
    logger.addHandler(file_handler)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(log_format)
    logger.addHandler(console_handler)
    logger.setLevel(logging_level)

    if inps.work_dir is None:
        work_dir = os.getcwd()
    else:
        work_dir = inps.work_dir
        if not os.path.exists(path=work_dir):
            logger.info(msg='create output folder: ' + work_dir)
            os.mkdir(path=work_dir)
    logger.info(msg='working directory: {}'.format(work_dir))

    downloadOSM(
        railway=inps.railway,
        highway=inps.highway,
        bridge=inps.bridge,
        work_dir=work_dir,
        out_file_name=inps.out_file_name,
        logger=logger,
        geom_file=inps.geom_file
    )


if __name__ == '__main__':
    main()
