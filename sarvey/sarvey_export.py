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

"""Console script for exporting data from SARvey format to GIS formats."""
import argparse
import logging
from logging import Logger
import sys
import time
import warnings
import os
from os.path import join, dirname, basename
import numpy as np
import pandas as pd
import geopandas as gpd
from pyproj import CRS, Transformer
from shapely import Point
from shapely.errors import ShapelyDeprecationWarning

from sarvey.config import loadConfiguration
from sarvey.console import showLogoSARvey
from sarvey.objects import Points, CoordinatesMap
import sarvey.utils as ut
from sarvey.geolocation import calculateGeolocationCorrection


warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)


def exportDataToGisFormat(*, file_path: str, output_path: str, input_path: str,
                          correct_geolocation: bool = False, no_timeseries: bool = False, logger: Logger):
    """Export data to GIS format (shp or gpkg).

    Parameters
    ----------
    file_path: str
        Path to the input file.
    output_path: str
        Path for writing output file.
    input_path: str
        Path to slcStack.h5 and geometryRadar.h5.
    correct_geolocation: bool
        Correct geolocation or not
    no_timeseries: bool
        Export time series data or not
    logger: Logger
        Logger handle.
    """
    point_obj = Points(file_path=file_path, logger=logger)

    point_obj.open(input_path=input_path)

    # todo: add corrected height to output
    # todo: add option to mask the output to e.g. linear infrastructures or other AOI

    vel, demerr, _, coherence, omega, _ = ut.estimateParameters(obj=point_obj, ifg_space=False)

    stc = ut.spatiotemporalConsistency(coord_map=point_obj.coord_map, phase=point_obj.phase,
                                       wavelength=point_obj.wavelength)

    point_obj.phase *= point_obj.wavelength / (4 * np.pi)  # in [m]

    # extract displacement
    defo_ts = np.zeros_like(point_obj.phase, dtype=np.float32)
    for i in range(point_obj.num_points):
        phase_topo = (point_obj.ifg_net_obj.pbase / (point_obj.slant_range[i] * np.sin(point_obj.loc_inc[i])) *
                      demerr[i])
        defo_ts[i, :] = point_obj.phase[i, :] - phase_topo

    # transform into meters
    defo_ts *= 1000  # in [mm]

    dates = ["D{}".format(date).replace("-", "") for date in point_obj.ifg_net_obj.dates]

    dates = dates[:point_obj.phase.shape[1]]  # remove dates which were not processed

    if no_timeseries:
        df_points = pd.DataFrame({})
    else:
        df_points = pd.DataFrame({date: () for date in dates})

    if correct_geolocation:
        logger.info("Calculate geolocation correction.")
        coord_correction = calculateGeolocationCorrection(path_geom=input_path,
                                                          point_obj=point_obj,
                                                          demerr=demerr,
                                                          logger=logger)
        coord_correction_norm = np.linalg.norm(coord_correction, axis=1)
        max_error_index = np.argmax(coord_correction_norm)
        logger.info(f"Maximum geolocation correction: {coord_correction_norm[max_error_index]:.1f} m "
                    f"corresponding to {demerr[max_error_index]:.1f} m DEM correction")
    else:
        coord_correction = 0
        logger.info("geolocation correction skipped.")

    coord_map = point_obj.coord_map
    coord_map += coord_correction

    # reconstruct Transverse Mercator projection
    coord_map_obj = CoordinatesMap(file_path=join(dirname(file_path), "coordinates_map.h5"), logger=logger)
    coord_map_obj.open()
    lon_0 = coord_map_obj.lon_0
    lat_0 = coord_map_obj.lat_0

    map_crs = CRS.from_dict({
        "proj": "tmerc",
        "lat_0": lat_0,
        "lon_0": lon_0,
        "k": 1,
        "x_0": 0,
        "y_0": 0,
        "datum": "WGS84",
        "units": "m",
        "no_defs": True
    })

    # construct Transverse Mercator to WGS84 transformer
    wgs_epsg = "4326"
    transformer = Transformer.from_crs(map_crs, wgs_epsg)

    lat, lon = transformer.transform(coord_map[:, 0], coord_map[:, 1])
    logger.debug(f"WGS84 Lon range: {np.min(lon):.6f} - {np.max(lon):.6f}")
    logger.debug(f"WGS84 Lat range: {np.min(lat):.6f} - {np.max(lat):.6f}")
    lonlat = np.vstack((lon, lat)).T

    logger.info('Construct output dataframe.')
    df_points['coord'] = (lonlat).tolist()
    df_points['coord'] = df_points['coord'].apply(Point)
    df_points.insert(0, 'point_id', point_obj.point_id.tolist())
    df_points.insert(1, 'velocity', vel * 1000)  # in [mm]
    df_points.insert(2, 'coherence', coherence)
    df_points.insert(3, 'omega', omega)
    df_points.insert(4, 'st_consistency', stc * 1000)  # in [mm]
    df_points.insert(5, 'dem_error', demerr)  # in [m]
    df_points.insert(6, 'dem', point_obj.height)  # in [m]

    df_points.columns = [col[:10] for col in df_points.columns]

    if not no_timeseries:
        for i, date in enumerate(dates):
            df_points[date] = defo_ts[:, i]

    gdf_points = gpd.GeoDataFrame(df_points, geometry='coord')
    gdf_points = gdf_points.set_crs(CRS.from_epsg(wgs_epsg))
    logger.info(msg="write to file.")
    gdf_points.to_file(output_path)


def createParser():
    """Create_parser."""
    parser = argparse.ArgumentParser(
        description="Export InSAR time series results from '.h5' to GIS data formats.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""EXAMPLE:
        sarvey_export outputs/p2_coh50_ts.h5 -o outputs/shp/p2_coh50.shp             # export time series to shapefile
        sarvey_export outputs/p2_coh50_ts.h5 -o outputs/shp/p2_coh50.gpkg            # export time series to geopackage
        sarvey_export outputs/p2_coh90_ts.h5 -o outputs/shp/p2_coh90.shp -g          # apply geolocation correction
        sarvey_export outputs/p2_coh90_ts.h5 -o outputs/shp/p2_coh90.shp -g -t       # skip time series data
        """)

    parser.add_argument('file_path', type=str, help='Path to input file.')

    parser.add_argument("-o", "--output_path", type=str, dest="output_path", default="",
                        help="Path to output file. If empty, the name of the input file will be used.")

    # parser.add_argument("-f", "--format", type=str, required=False, metavar="FILE", dest="format",
    #                     help="Output file format (if not already specified within '-o'). Can be 'shp', 'gpkg',
    #                     'csv'.")

    parser.add_argument("-l", "--log_dir", type=str, required=False, metavar="FILE", dest="log_dir",
                        default="logfiles/", help="Logfile directory (default: 'logfiles/')")

    parser.add_argument('-w', '--workdir', default=None, dest="workdir",
                        help='Working directory (default: current directory).')

    parser.add_argument('-g', '--correct_geo', default=False, action="store_true", dest="correct_geolocation",
                        help='Correct Geolocation (default: False).')

    parser.add_argument('-t', '--no-time-series', default=False, action="store_true", dest="no_timeseries",
                        help='Do not export time series (default: False).')

    return parser


def main(iargs=None):
    """Run Main.

    :param iargs:
    """
    parser = createParser()
    args = parser.parse_args(iargs)

    log_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(log_format)
    logger.addHandler(console_handler)
    logging_level = logging.getLevelName("INFO")
    logger.setLevel(logging_level)

    if args.workdir is None:
        args.workdir = os.path.abspath(os.path.curdir)
    else:
        logger.info(msg="Working directory: {}".format(args.workdir))

    args.log_dir = join(args.workdir, args.log_dir)
    current_datetime = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
    log_filename = f"sarvey_export_{current_datetime}.log"

    if not os.path.exists(args.log_dir):
        os.mkdir(args.log_dir)
    file_handler = logging.FileHandler(filename=join(args.log_dir, log_filename))
    file_handler.setFormatter(log_format)
    logger.addHandler(file_handler)

    showLogoSARvey(logger=logger, step="Export results")

    # read config file to retrieve location of inputs
    config_file_path = os.path.abspath(join(args.workdir, dirname(args.file_path), "config.json"))

    if not os.path.exists(config_file_path):
        # check if any config file is available in upper directory (backward compatibility)
        files = np.array([os.path.abspath(f) for f in os.listdir(join(dirname(config_file_path), ".."))
                          if os.path.isfile(f)])
        potential_configs = np.array([(basename(f).split(".")[-1] == "json") and ("config" in basename(f))
                                      for f in files])
        if potential_configs[potential_configs].shape[0] == 0:
            raise FileNotFoundError(f"Backup configuration file not found: {config_file_path}!")
        else:
            logger.warning(msg=f"Backup configuration file not found: {config_file_path}!")
            logger.warning(msg=f"Other configuration files automatically detected: {files[potential_configs]}!")
            logger.warning(msg=f"Automatically selected configuration file: {files[potential_configs][0]}!")
            config_file_path = files[potential_configs][0]

    config = loadConfiguration(path=config_file_path)

    # create output directory
    if args.output_path == "":
        output_dir = args.workdir
        output_fname = basename(args.file_path).split(".")[-2]
        output_format = "shp"
        args.output_path = join(output_dir, output_fname + "." + output_format)
    else:
        output_dir = join(args.workdir, dirname(args.output_path))
        output_fname = basename(args.output_path)
        name_splitted = output_fname.split(".")
        if len(name_splitted) == 1:
            args.output_path = join(output_dir, output_fname + ".shp")  # use shapefile as default format
        elif len(name_splitted) == 2:
            output_format = name_splitted[-1]  # use specified format
            if (output_format != "shp") and (output_format != "gpkg"):
                logger.error(msg=f"Output format not supported: {output_format}!")
                raise ValueError
            logger.info(msg=f"Detected output format: {output_format}.")
            args.output_path = join(output_dir, output_fname)
        else:
            logger.error(msg=f"Output format was not recognized! {output_fname}")
            raise ValueError

    logger.info(msg=f"Output file: {args.output_path}")

    # specify geolocation status
    logger.info(msg=f"Correct geolocation: {args.correct_geolocation}")

    # specify time series flag
    logger.info(msg=f"Export time series data: {not args.no_timeseries}")

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    exportDataToGisFormat(file_path=args.file_path, output_path=args.output_path,
                          input_path=config.general.input_path,
                          correct_geolocation=args.correct_geolocation, no_timeseries=args.no_timeseries,
                          logger=logger)


if __name__ == '__main__':
    main()
