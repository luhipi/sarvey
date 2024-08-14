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

"""Generate mask from shape file."""
import argparse
import os
import sys
import time
from os.path import join
import PIL.Image as Image
import PIL.ImageDraw as ImageDraw
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy import spatial
import logging
from logging import Logger
import geopandas as gpd

from mintpy.utils import writefile, ptime, utils

from sarvey.osm_utils import getSpatialExtend

try:
    matplotlib.use('TkAgg')
except ImportError as e:
    print(e)

EXAMPLE = """Example:
  sarvey_mask path/to/file.shp --geom ./geometryRadar.h5 --width 6 -o mask_infra.h5
"""


def create_parser():
    """Create_parser."""
    parser = argparse.ArgumentParser(
        description='Create transport infrastructure mask from shp-file.',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=EXAMPLE)

    parser.add_argument(dest='input_file', help='path to input shp-file.')

    parser.add_argument('-w', '--work_dir', dest='work_dir', default=None,
                        help='absolute path to working directory\n' +
                             '(default: current directory).')

    parser.add_argument('--geom', dest='geom_file', default=None,
                        help='path to existing geometryRadar.h5 file.')

    parser.add_argument('--width', dest='width', default=6, type=int,
                        help='Width of the mask in pixel (default: 6).')

    parser.add_argument('-o', dest='out_file_name', default='mask_infra.h5',
                        help="name of output file. (default: 'mask_infra.h5').")

    return parser


class Node:
    """Define simple class for a node at a road (similar to overpy.Node)."""

    def __init__(self, *, lat: float = None, lon: float = None):
        """Init."""
        self.lat = lat
        self.lon = lon


class CoordinateSearch:
    """CoordinateSearch."""

    def __init__(self):
        """Init."""
        self.search_tree = None
        self.yidx = None
        self.xidx = None
        self.lon = None
        self.lat = None
        self.coord = None

    def createSearchTree(self, *, coord: utils.coordinate, logger: Logger):
        """Create search tree.

        Parameters
        ----------
        coord: utils.coordinate
            Coordinates
        logger: Logger
            Logging handler.
        """
        self.coord = coord
        logger.info('create kd-tree for efficient search...')

        if self.coord.lut_y is None or self.coord.lut_x is None:
            self.coord.open()
            lat, lon = self.coord.read_lookup_table(print_msg=False)
            self.lat = lat.ravel()
            self.lon = lon.ravel()

        # create the 2D coordinate arrays for fast indexing
        x = np.arange(self.coord.lut_x.shape[1])
        y = np.arange(self.coord.lut_x.shape[0])
        xx, yy = np.meshgrid(x, y)
        self.xidx = xx.ravel()
        self.yidx = yy.ravel()

        start_time = time.time()
        self.search_tree = spatial.KDTree(data=np.array([self.lon, self.lat]).transpose())

        logger.info('... done.')
        m, s = divmod(time.time() - start_time, 60)
        logger.debug(f"time used: {m:02.0f} mins {s:02.1f} secs.")

    def getMeanDistanceBetweenPixels(self):
        """Compute mean distance between adjacent pixels."""
        distances = self.search_tree.query([self.lon[0], self.lat[0]], k=10)[0]
        mean_dist = np.mean(distances[1:])
        return mean_dist

    def getNearestNeighbour(self, *, node: Node):
        """Query the kd-tree for the nearest neighbour.

        :param node: Node object
        """
        # find nearest neighbour
        dist, idx = self.search_tree.query([node.lon, node.lat])
        found_node = Node(lat=self.lat[idx], lon=self.lon[idx])
        # return index of NN in radar coordinates
        return dist, (self.yidx[idx], self.xidx[idx]), found_node


def findLastRoadPixel(*, csearch: CoordinateSearch, cur_node: Node, prev_node: Node, dist_thrsh: float):
    """Find the index of the last road pixel that is within the image extend.

    Idea: the pixel with the shortest distance to the current node of a road is not necessarily on the road, if the
    current node is outside the image extend. Split the road in further linear parts and find the last road pixel
    recursively that is still inside the image.
    Hint: all nodes are instances from class Node

    Parameters
    ----------
    csearch: CoordinateSearch
        Search tree for efficient spatial search of the coordinate of a pixel in the radar image.
    cur_node: Node
        Current node of the road that is outside the image extend.
    prev_node: Node
        Previous node of the road that is inside the image extend.
    dist_thrsh: float
        Distance threshold for stop criterion (derived from average distance between two pixels in the image).

    Returns
    -------
    node_idx: int
        Node of the pixel which is the last pixel on the road inside the image.
    """
    # create a new node at half of the road distance between previous and current node
    mid_lat = cur_node.lat + (cur_node.lat - prev_node.lat) / 2
    mid_lon = cur_node.lon + (cur_node.lon - prev_node.lon) / 2
    mid_node = Node(lat=mid_lat, lon=mid_lon)

    dist, node_idx = csearch.getNearestNeighbour(node=mid_node)[0:2]
    if dist < dist_thrsh:
        return node_idx
    else:
        node_idx = findLastRoadPixel(csearch=csearch, cur_node=cur_node, prev_node=prev_node, dist_thrsh=dist_thrsh)
    return node_idx


def euclDist(*, node1: Node, node2: Node):
    """Compute the euclidean distance between two nodes."""
    return np.sqrt((node1.lat - node2.lat) ** 2 + (node1.lon - node2.lon) ** 2)


def computeLastRoadPixel(*, cur_node: Node, prev_node: Node, found_node: Node):
    """Compute the location of the pixel at the border of the radar image that is part of the road.

    Parameters
    ----------
    cur_node: Node
        Current node of the road.
    prev_node: Node
        Previous node of the road.
    found_node: Node
        Found node of the road.

    Returns
    -------
    new_lon: float
        Longitude of the pixel at the border of the radar image that is part of the road.
    new_lat: float
        Latitude of the pixel at the border of the radar image that is part of the road.
    """
    a = euclDist(node1=prev_node, node2=found_node)
    b = euclDist(node1=cur_node, node2=found_node)
    c = euclDist(node1=prev_node, node2=cur_node)
    alpha = np.arccos((- a ** 2 + b ** 2 + c ** 2) / (2 * b * c))
    d = b / np.sin(np.pi / 2 - alpha)
    new_lat = cur_node.lat + (prev_node.lat - cur_node.lat) / c * d
    new_lon = cur_node.lon + (prev_node.lon - cur_node.lon) / c * d
    return new_lon, new_lat


def convertToRadarCoordPolygon(*, gdf_infra: gpd.geodataframe, csearch: CoordinateSearch, logger: Logger):
    """Convert Polygon to a mask in shape of radar image.

    Parameters
    ----------
    gdf_infra: gpd.geodataframe
        The queried infrastructures containing polygons.
    csearch: CoordinateSearch
        The coordinate search object.
    logger: Logger
        Logging handler.

    Returns
    -------
    img_np: np.ndarray
        Mask image.
    """
    # create a new image
    logger.info('create mask image...')
    img_pil = Image.new(mode="1",
                        size=(int(csearch.coord.src_metadata['LENGTH']), int(csearch.coord.src_metadata['WIDTH'])))
    img_pil_draw = ImageDraw.Draw(im=img_pil)

    num_ways = gdf_infra.shape[0]
    way_iter = 0
    prog_bar = ptime.progressBar(maxValue=num_ways)

    dist_thrsh = 1.3 * csearch.getMeanDistanceBetweenPixels()
    lines = [geom.boundary.coords for geom in gdf_infra.geometry if geom is not None]
    way_list = list()
    for coo in lines:
        way_list.append([Node(lat=point[1], lon=point[0]) for point in coo])

    # plt.ion()
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel("lon")
    ax.set_ylabel("lat")
    lat, lon = csearch.coord.read_lookup_table(print_msg=False)
    ax.plot([lon[0, 0], lon[-1, 0]], [lat[0, 0], lat[-1, 0]], '-k')
    ax.plot([lon[0, 0], lon[0, -1]], [lat[0, 0], lat[0, -1]], '-k')
    ax.plot([lon[0, -1], lon[-1, -1]], [lat[0, -1], lat[-1, -1]], '-k')
    ax.plot([lon[-1, 0], lon[-1, -1]], [lat[-1, 0], lat[-1, -1]], '-k')
    # ax.plot(lon.ravel(), lat.ravel(), '.k', markersize=0.5)

    while way_iter < num_ways:
        way = way_list[way_iter]
        poly_line_way = []

        # perform a preliminary search to check if polygon is partly outside image extend
        outside = np.zeros(len(way))
        for i in range(len(way)):
            cur_node = way[i]

            # convert node coordinates (lat, lon) to image coordinates
            dist, _, _ = csearch.getNearestNeighbour(node=cur_node)

            # check if node is outside the image
            if dist > dist_thrsh:
                outside[i] = 1

        if np.sum(outside) == 0:  # all road nodes inside image extend
            for i in range(len(way)):
                cur_node = way[i]

                # convert node coordinates (lat, lon) to image coordinates
                dist, node_idx, found_node = csearch.getNearestNeighbour(node=cur_node)

                # Fill list of current way with node coordinates
                poly_line_way.append(node_idx)
                ax.plot(cur_node.lon, cur_node.lat, '*k')
                ax.plot(found_node.lon, found_node.lat, 'ok')

        else:  # some polygon nodes outside image extend
            if np.sum(outside) == outside.size:  # all nodes outside, skip
                way_iter += 1
                continue

            # polygon nodes partly inside and partly outside
            prev_p = outside[-2] == 1  # last point == first point (due to closed polygon). Select second last.
            for i in range(outside.shape[0]):
                cur_p = outside[i] == 1
                cur_node = way[i]

                # convert node coordinates (lat, lon) to image coordinates
                dist, node_idx, found_node = csearch.getNearestNeighbour(node=cur_node)

                # check if transition happens
                #   yes: check if current point is inside or outside
                #        if outside: find transition point, but do not add current point
                #        if inside: find transition point, then add current point
                #   no:  if point inside: add point
                #        if point outside: skip point

                if not (prev_p == cur_p):  # transition
                    stored_idx = None
                    if i - 1 < 0:
                        prev_node = way[-2]
                    else:
                        prev_node = way[i - 1]

                    if cur_p:  # transition: in -> out
                        # find transition point, but do not add current point
                        ax.plot(cur_node.lon, cur_node.lat, '*y')

                    if prev_p:  # transition: out -> in
                        # find and add transition point, then add current point.
                        stored_idx = node_idx  # store current point for adding it later.
                        ax.plot(cur_node.lon, cur_node.lat, '*r')  # plot now, because variables will be overwritten
                        # the 'found_node' has to be computed from the last point outside, i.e. from 'prev_node'
                        ax.plot(found_node.lon, found_node.lat, 'or')
                        _, _, found_node = csearch.getNearestNeighbour(node=prev_node)

                    new_lon, new_lat = computeLastRoadPixel(
                        cur_node=cur_node,
                        prev_node=prev_node,
                        found_node=found_node
                    )

                    dist, node_idx, found_node = csearch.getNearestNeighbour(node=Node(lon=new_lon, lat=new_lat))
                    ax.plot(cur_node.lon, cur_node.lat, '*b')
                    ax.plot(found_node.lon, found_node.lat, 'ob')
                    ax.plot(new_lon, new_lat, '+b')

                    # add the transition point
                    poly_line_way.append(node_idx)
                    if prev_p:  # transition: out -> in
                        # transition point found and added, now add stored current point.
                        poly_line_way.append(stored_idx)
                    prev_p = cur_p  # prepare for next iteration

                elif cur_p:  # no transition, current point is outside -> do not add point
                    ax.plot(cur_node.lon, cur_node.lat, '*y')
                    prev_p = cur_p  # prepare for next iteration

                else:  # no transition, current points is inside -> add point
                    ax.plot(cur_node.lon, cur_node.lat, '*r')
                    ax.plot(found_node.lon, found_node.lat, 'or')
                    poly_line_way.append(node_idx)
                    prev_p = cur_p  # prepare for next iteration

        prog_bar.update(value=way_iter + 1, every=10, suffix='{}/{} polygons'.format(way_iter + 1, num_ways))

        # if first point is outside image, the polygon will not be closed. However, it still works to create a polygon.
        img_pil_draw.polygon(poly_line_way, fill=255)
        # plt.figure()
        # plt.imshow(np.array(img_pil.getdata()).reshape(img_pil.size[1], img_pil.size[0]).astype(int))

        way_iter += 1

    img_np = np.array(img_pil.getdata()).reshape(img_pil.size[1], img_pil.size[0]).astype(int)
    return img_np


def convertToRadarCoord(*, gdf_infra: gpd.geodataframe, csearch: CoordinateSearch, width: int, logger: Logger):
    """Convert Polyline to a mask in shape of radar image. Apply a buffer of size 'width' in pixels.

    Parameters
    ----------
    gdf_infra: gpd.geodataframe
        The queried infrastructures containing polygons.
    csearch: CoordinateSearch
        The coordinate search object.
    width: int
        Width of the mask in pixel.
    logger: Logger
        Logging handler.

    Returns
    -------
    img_np: np.ndarray
        Mask image.
    """
    # create a new image
    logger.info('create mask image...')
    img_pil = Image.new(mode="1",
                        size=(int(csearch.coord.src_metadata['LENGTH']), int(csearch.coord.src_metadata['WIDTH'])))
    img_pil_draw = ImageDraw.Draw(im=img_pil)

    num_roads = gdf_infra.shape[0]
    prog_bar = ptime.progressBar(maxValue=num_roads)

    dist_thrsh = 1.3 * csearch.getMeanDistanceBetweenPixels()
    lines = [ls.coords for ls in gdf_infra.geometry if ls is not None]  # enables to append to list
    way_list = list()
    for coo in lines:
        way_list.append([Node(lat=point[1], lon=point[0]) for point in coo])

    num_ways = len(way_list)  # changes during iteration
    way_iter = 0

    # plt.ion()
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel("lon")
    ax.set_ylabel("lat")
    lat, lon = csearch.coord.read_lookup_table(print_msg=False)
    ax.plot([lon[0, 0], lon[-1, 0]], [lat[0, 0], lat[-1, 0]], '-k')
    ax.plot([lon[0, 0], lon[0, -1]], [lat[0, 0], lat[0, -1]], '-k')
    ax.plot([lon[0, -1], lon[-1, -1]], [lat[0, -1], lat[-1, -1]], '-k')
    ax.plot([lon[-1, 0], lon[-1, -1]], [lat[-1, 0], lat[-1, -1]], '-k')
    # ax.plot(lon.ravel(), lat.ravel(), '.k', markersize=0.5)

    while way_iter < num_ways:
        way = way_list[way_iter]
        poly_line_way = []

        # perform a preliminary search to check if road is partly outside image extend
        outside = np.zeros(len(way))
        for i in range(len(way)):
            cur_node = way[i]

            # convert node coordinates (lat, lon) to image coordinates
            dist, _, _ = csearch.getNearestNeighbour(node=cur_node)

            # check if node is outside the image
            if dist > dist_thrsh:
                outside[i] = 1

        if np.sum(outside) == 0:  # all road nodes inside image extend
            for i in range(len(way)):
                cur_node = way[i]

                # convert node coordinates (lat, lon) to image coordinates
                dist, node_idx, found_node = csearch.getNearestNeighbour(node=cur_node)

                # Fill list of current way with node coordinates
                poly_line_way.append(node_idx)
                ax.plot(cur_node.lon, cur_node.lat, '*k')
                ax.plot(found_node.lon, found_node.lat, 'ok')

        else:  # some road nodes outside image extend
            if np.sum(outside) == outside.size:  # all nodes outside, skip
                way_iter += 1
                continue
            # split the way into sub parts based on in-out / out-in transition
            # find first node inside the image
            first_inside_idx = np.where(outside == 0)[0][0]
            if first_inside_idx > 0:  # this is a transition into the image
                start_idx = first_inside_idx - 1
            else:
                start_idx = first_inside_idx

            # find first node which is again outside the image
            outside_idx = np.where(outside[first_inside_idx:] == 1)[0]
            if outside_idx.size == 0:  # no more transition to outside the image
                stop_idx = len(way)
            else:
                stop_idx = outside_idx[0] + first_inside_idx + 1
                if stop_idx != len(way):  # split the current way and add a new way at the end of the way_list
                    # to handle it later
                    way_list.append(way[stop_idx:])
                    num_ways += 1

            for i in range(start_idx, stop_idx):
                cur_node = way[i]

                # convert node coordinates (lat, lon) to image coordinates
                dist, node_idx, found_node = csearch.getNearestNeighbour(node=cur_node)

                if dist > dist_thrsh:
                    if i == start_idx:  # there is no previous node, but a next node.
                        prev_node = way[i + 1]
                    else:
                        prev_node = way[i - 1]
                    new_lon, new_lat = computeLastRoadPixel(cur_node=cur_node, prev_node=prev_node,
                                                            found_node=found_node)
                    dist, node_idx, found_node = csearch.getNearestNeighbour(node=Node(lon=new_lon, lat=new_lat))
                    ax.plot(cur_node.lon, cur_node.lat, '*b')
                    ax.plot(found_node.lon, found_node.lat, 'ob')
                    ax.plot(new_lon, new_lat, '+b')
                else:
                    ax.plot(cur_node.lon, cur_node.lat, '*r')
                    ax.plot(found_node.lon, found_node.lat, 'or')
                # Fill list of current way with node coordinates
                poly_line_way.append(node_idx)

        prog_bar.update(value=way_iter + 1, every=10, suffix='{}/{} road segments'.format(way_iter + 1, num_roads))

        img_pil_draw.line(poly_line_way, fill=255, width=width)
        # img_pil_draw.polygon(poly_line_way, fill=255)

        way_iter += 1

    img_np = np.array(img_pil.getdata()).reshape(img_pil.size[1], img_pil.size[0]).astype(int)
    return img_np


def saveMask(*, work_dir: str, mask: np.ndarray, atr: dict, out_file_name: str):
    """Save the mask to 'maskRoads.h5'.

    Parameters
    ----------
    work_dir: str
        Working directory.
    mask: np.ndarray
        Mask image.
    atr: dict
        Metadata data, e.g. from the geometryRadar.h5 file.
    out_file_name: str
        Output file name.
    """
    # create the right attributes
    ds_dict = dict()
    ds_dict['mask'] = mask.transpose().astype('float32')
    atr["FILE_TYPE"] = "mask"

    writefile.write(datasetDict=ds_dict, out_file=os.path.join(work_dir, out_file_name), metadata=atr)


def createMask(*, input_file: str, width: int, work_dir: str, out_file_name: str, geom_file: str,
               logger: logging.Logger):
    """Create a mask for the radar image from a shapefile containing lines or polygons.

    Parameters
    ----------
    input_file: str
        Path to input file.
    width: int
        Width of the mask in pixel. Applied to the lines only.
    work_dir: str
        Working directory.
    out_file_name: str
        Output file name.
    geom_file: str
        Path to geometryRadar.h5 file.
    logger: logging.Logger
        Logging handler.
    """
    logger.info("Start creating mask file based on openstreetmap data.")

    # get bounding box
    _, _, _, coord, atr = getSpatialExtend(geom_file=geom_file, logger=logger)

    # create search tree
    csearch = CoordinateSearch()
    csearch.createSearchTree(coord=coord, logger=logger)

    logger.info(f"Read from input file: {input_file}")
    gdf_infra = gpd.read_file(input_file)

    if gdf_infra.geometry[0].geom_type == "LineString":
        mask_img = convertToRadarCoord(gdf_infra=gdf_infra, csearch=csearch, width=width, logger=logger)

    elif gdf_infra.geometry[0].geom_type == "Polygon":
        mask_img = convertToRadarCoordPolygon(gdf_infra=gdf_infra, csearch=csearch, width=width, logger=logger)
    else:
        logger.error(f"Geometry type is {gdf_infra.geometry[0].geom_type}."
                         f"Only 'LineString' and 'Polygon' supported!")
        raise TypeError

    if '.h5' not in out_file_name:
        out_file_name += ".h5"
    saveMask(work_dir=work_dir, mask=mask_img, atr=atr, out_file_name=out_file_name)

    logger.info("Masking finished.")


def main(iargs=None):
    """Create mask from lines or polygons given in geographic coordinates (EPSG:4326). Input as shp or gpkg."""
    # check input
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # initiate logger
    logging_level = logging.getLevelName('DEBUG')

    log_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)

    current_datetime = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
    log_filename = f"sarvey_mask_{current_datetime}.log"
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
            logger.info('create output folder: ' + work_dir)
            os.mkdir(path=work_dir)
    logger.info(f"working directory: {work_dir}")

    input_file = join(work_dir, inps.input_file)
    out_file_name = join(work_dir, inps.out_file_name)

    createMask(
        input_file=input_file,
        width=inps.width,
        work_dir=work_dir,
        out_file_name=out_file_name,
        logger=logger,
        geom_file=inps.geom_file
    )


if __name__ == '__main__':
    main()
