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

"""Objects module for SARvey."""
import os
from os.path import join, dirname, exists, basename
from typing import Optional
import h5py
import matplotlib.pyplot as plt
import numpy as np
from pyproj import Proj, CRS
from pyproj.aoi import AreaOfInterest
from pyproj.database import query_utm_crs_info
from logging import Logger

from miaplpy.objects.slcStack import slcStack
from mintpy.utils import readfile
from mintpy.utils.plot import auto_flip_direction

from sarvey.ifg_network import IfgNetwork


class AmplitudeImage:
    """AmplitudeImage."""

    def __init__(self, *, file_path: str):
        """Init.

        Parameters
        ----------
        file_path: str
            path to filename
        """
        self.width = None
        self.length = None
        self.file_path = file_path
        self.background_map = None
        self.orbit_direction = None

    def prepare(self, *, slc_stack_obj: slcStack, img: np.ndarray, logger: Logger):
        """Read the SLC stack, compute the mean amplitude image and store it into a file.

        Parameters
        ----------
        slc_stack_obj: slcStack
            object of class slcStack from MiaplPy
        img: np.ndarray
            amplitude image, e.g. the mean over time
        logger: Logger
            Logging handler
        """
        self.orbit_direction = slc_stack_obj.metadata["ORBIT_DIRECTION"]
        self.length = slc_stack_obj.length
        self.width = slc_stack_obj.width

        self.background_map = img

        logger.info(msg="write data to {}...".format(self.file_path))

        if exists(self.file_path):
            os.remove(self.file_path)

        with h5py.File(self.file_path, 'w') as f:
            f.create_dataset('background_map', data=self.background_map)
            f.attrs["ORBIT_DIRECTION"] = self.orbit_direction
            f.attrs["LENGTH"] = self.length
            f.attrs["WIDTH"] = self.width

    def open(self):
        """Open."""
        # print("read from {}".format(self.file_path))

        with h5py.File(self.file_path, 'r') as f:
            self.background_map = f["background_map"][:]
            self.orbit_direction = f.attrs["ORBIT_DIRECTION"]
            self.length = f.attrs["LENGTH"]
            self.width = f.attrs["WIDTH"]

    def plot(self, *, ax: plt.Axes = None, logger: Logger):
        """Plot the mean amplitude image as a background map.

        Parameters
        ----------
        ax: plt.Axes
            axes for plotting (default: None, a new figure will be created).
        logger: Logger
            Logging handler.

        Return
        ------
        ax: plt.Axes
            axes object.
        """
        if self.background_map is None:
            try:
                self.open()
            except OSError as e:
                logger.error(msg="Could not open file: {}".format(e))
                fig = plt.figure(figsize=(15, 5))
                ax = fig.add_subplot()
                logger.error(msg="Orbit direction not available.")
                return ax

        if ax is None:
            fig = plt.figure(figsize=(15, 5))
            ax = fig.add_subplot()
        ax.imshow(self.background_map, cmap=plt.cm.get_cmap("gray"))
        meta = {"ORBIT_DIRECTION": self.orbit_direction}
        auto_flip_direction(meta, ax=ax, print_msg=False)

        ax.set_xlabel("Range")
        ax.set_ylabel("Azimuth")

        return ax


class CoordinatesUTM:
    """Coordinates in UTM for all pixels in the radar image."""

    def __init__(self, *, file_path: str, logger: Logger):
        """Init.

        Parameters
        ----------
        file_path: str
            path to filename
        logger: Logger
            Logging handler.
        """
        self.file_path = file_path
        self.coord_utm = None
        self.logger = logger

    def prepare(self, *, input_path: str):
        """Read the slc stack, computes the mean amplitude image and stores it into a file.

        Parameters
        ----------
        input_path: str
            path to slcStack.h5 file.
        """
        log = self.logger
        lat = readfile.read(input_path, datasetName='latitude')[0]
        lon = readfile.read(input_path, datasetName='longitude')[0]

        log.info(msg="Transform coordinates from latitude and longitude (WGS84) to North and East (UTM).")
        # noinspection PyTypeChecker
        utm_crs_list = query_utm_crs_info(
            datum_name="WGS 84",
            area_of_interest=AreaOfInterest(
                west_lon_degree=np.nanmin(lon.ravel()),
                south_lat_degree=np.nanmin(lat.ravel()),
                east_lon_degree=np.nanmax(lon.ravel()),
                north_lat_degree=np.nanmax(lat.ravel())),
            contains=True)
        utm_crs = CRS.from_epsg(utm_crs_list[0].code)
        lola2utm = Proj(utm_crs)
        self.coord_utm = np.array(lola2utm(lon, lat))

        log.info(msg="write data to {}...".format(self.file_path))

        if exists(self.file_path):
            os.remove(self.file_path)

        with h5py.File(self.file_path, 'w') as f:
            f.create_dataset('coord_utm', data=self.coord_utm)

    def open(self):
        """Open."""
        with h5py.File(self.file_path, 'r') as f:
            self.coord_utm = f["coord_utm"][:]


class BaseStack:
    """Class for 3D image-like data stacks."""

    def __init__(self, *, file: str = None, logger: Logger):
        """Init.

        Parameters
        ----------
        file: str
            path to filename
        logger: Logger
            Logging handler.
        """
        self.file = file
        self.logger = logger
        self.metadata = None
        self.num_time = None
        self.length = None
        self.width = None
        self.f = None

    def close(self, *, print_msg: bool = True):
        """Close."""
        try:
            self.f.close()
            if print_msg:
                self.logger.info(msg='close file: {}'.format(basename(self.file)))
        except Exception as e:
            self.logger.exception(msg=e)
            pass
        return None

    def getShape(self, *, dataset_name: str):
        """Open file and read shape of dataset."""
        with h5py.File(self.file, 'r') as f:
            dshape = f[dataset_name].shape
        return dshape

    def read(self, *, dataset_name: str, box: Optional[tuple] = None, print_msg: bool = True):
        """Read dataset from slc file.

        Parameters
        ----------
        dataset_name: str
            name of dataset
        box: tuple
            tuple of 4 int, indicating x0,y0,x1,y1 of range, or
            tuple of 6 int, indicating x0,y0,z0,x1,y1,z1 of range
        print_msg: bool
            print message.

        Returns
        -------
        data: np.ndarray
            2D or 3D dataset
        """
        if print_msg:
            self.logger.info(msg='reading box {} from file: {} ...'.format(box, self.file))

        with h5py.File(self.file, 'r') as f:
            self.metadata = dict(f.attrs)

            ds = f[dataset_name]
            if len(ds.shape) == 3:
                self.length, self.width, self.num_time = ds.shape
            else:
                self.length, self.width = ds.shape

            # Get Index in space/2_3 dimension
            if box is None:
                box = [0, 0, self.width, self.length]

            if len(ds.shape) == 3:
                if len(box) == 4:
                    data = ds[box[1]:box[3], box[0]:box[2], :]
                if len(box) == 6:
                    data = ds[box[1]:box[4], box[0]:box[3], box[2]:box[5]]
            else:
                if len(box) == 6:
                    raise IndexError("Cannot read 3D box from 2D data.")
                data = ds[box[1]:box[3], box[0]:box[2]]

        for key, value in self.metadata.items():
            try:
                self.metadata[key] = value.decode('utf8')
            except Exception:
                self.metadata[key] = value
        return data

    def prepareDataset(self, dataset_name: str, dshape: tuple, dtype: object,
                       metadata: Optional[dict], mode: str = "w", chunks: [tuple, bool] = True):
        """PrepareDataset. Creates a dataset in file with specified size without writing any data.

        Parameters
        ----------
        dataset_name: str
            name of dataset.
        dshape: tuple
            shape of dataset.
        dtype: object
            data type of dataset.
        metadata: dict
            metadata of dataset (e.g. WAVELENGTH, ORBIT_DIRECTION, etc.). Usually the same as in slcStack.h5.
        mode: str
            open mode ('w' for writing new file or 'a' for appending to existing file).
        chunks: tuple
            chunk size ('True'/'False' or tuple specifying the dimension of the chunks)
        """
        with h5py.File(self.file, mode) as f:
            self.logger.info(msg="Prepare dataset: {d:<25} of {t:<25} in size of {s}".format(
                d=dataset_name,
                t=str(dtype),
                s=dshape))

            f.create_dataset(dataset_name,
                             shape=dshape,
                             dtype=dtype,
                             chunks=chunks)

            # write attributes
            metadata = dict(metadata)
            for key in metadata.keys():
                f.attrs[key] = metadata[key]

        return

    def writeToFileBlock(self, *, data: np.ndarray, dataset_name: str, block: Optional[tuple] = None, mode: str = 'a',
                         print_msg: bool = True):
        """Write data to existing HDF5 dataset in disk block by block.

        Parameters
        ----------
        data: np.ndarray
            1/2/3D matrix.
        dataset_name: str
            dataset name.
        block: list
            the list can contain 2, 4 or 6 integers indicating: [zStart, zEnd, yStart, yEnd, xStart, xEnd].
        mode: str
            open mode ('w' for writing new file or 'a' for appending to existing file).
        print_msg: bool
            print message.

        Returns
        --------
        file: str
            path to file
        """
        if block is None:
            # data shape
            if isinstance(data, list):
                shape = (len(data),)
            else:
                shape = data.shape

            if len(shape) == 1:
                block = [0, shape[0]]
            elif len(shape) == 2:
                block = [0, shape[0],
                         0, shape[1]]
            elif len(shape) == 3:
                block = [0, shape[0],
                         0, shape[1],
                         0, shape[2]]

        with h5py.File(self.file, mode) as f:

            if print_msg:
                self.logger.info(msg="writing dataset /{:<25} block: {}".format(dataset_name, block))
            if len(block) == 6:
                f[dataset_name][block[0]:block[1],
                                block[2]:block[3],
                                block[4]:block[5]] = data

            elif len(block) == 4:
                f[dataset_name][block[0]:block[1],
                                block[2]:block[3]] = data

            elif len(block) == 2:
                f[dataset_name][block[0]:block[1]] = data

        return self.file

    def writeToFile(self, *, data: np.ndarray, dataset_name: str, metadata: Optional[dict] = None, mode: str = 'a',
                    chunks: [tuple, bool] = True):
        """Write the whole dataset to the file (not block-by-block).

        Parameters
        ----------
        data: np.ndarray
            3D data array.
        dataset_name: str
            name of dataset.
        metadata: dict
            metadata of dataset (e.g. WAVELENGTH, ORBIT_DIRECTION, etc.). Usually the same as in slcStack.h5.
        mode: str
            mode for opening the h5 file (e.g. write: 'w' or append: 'a')
        chunks: tuple
            chunk size ('True'/'False' or tuple specifying the dimension of the chunks)
        """
        # 3D dataset
        self.logger.info(msg='create HDF5 file: {} with w mode'.format(self.file))
        self.f = h5py.File(self.file, mode)
        if dataset_name not in self.f:
            self.logger.info(msg='create dataset /{n} of {t:<10} in size of {s}.'.format(n=dataset_name,
                                                                                         t=str(data.dtype),
                                                                                         s=data.shape))
            self.f.create_dataset(dataset_name, data=data, chunks=chunks)
        else:
            self.logger.info(msg='overwrite dataset /{n} of {t:<10} in size of {s}.'.format(n=dataset_name,
                                                                                            t=str(data.dtype),
                                                                                            s=data.shape))
            self.f[dataset_name] = data

        # Attributes
        if metadata is not None:
            metadata = dict(metadata)
            for key, value in metadata.items():
                self.f.attrs[key] = str(value)

        self.f.close()
        self.logger.info(msg='finished writing to {}'.format(self.file))
        return


class Points:
    """Points class for storing information about the selected scatterers."""

    file_path: str
    point_id: np.array
    coord_xy: np.array
    num_points: int
    phase: np.array
    wavelength: float
    length: int
    width: int
    times: None

    # etc.

    def __init__(self, *, file_path: str, logger: Logger):
        """Init.

        Parameters
        ----------
        file_path: str
             ath to filename
        logger: Logger
            Logging handler.
        """
        self.ifg_net_obj = IfgNetwork()  # use parent class here which doesn't know and care about 'star' or 'sb'
        self.coord_utm = None
        self.coord_lalo = None
        self.height = None
        self.slant_range = None
        self.loc_inc = None
        self.file_path = file_path
        self.logger = logger

    def prepare(self, *, point_id: np.ndarray, coord_xy: np.ndarray, path_inputs: str):
        """Assign point_id and radar coordinates to the object.

        Store the point_id and radar coordinates of the scatterers in the object (not file) and read further
        attributes from external files (ifg_network.h5, slcStack.h5, geometryRadar.h5, coordinates_utm.h5).

        Parameters
        ----------
        point_id: np.ndarray
            point_id of the scatterers.
        coord_xy: np.ndarray
            radar coordinates of the scatterers.
        path_inputs: str
            path to input files (slcStack.h5, geometryRadar.h5).
        """
        self.point_id = point_id
        self.coord_xy = coord_xy
        self.num_points = self.coord_xy.shape[0]
        self.phase = None
        self.openExternalData(path_inputs=path_inputs)

    def writeToFile(self):
        """Write data to .h5 file (num_points, coord_xy, point_id, phase)."""
        self.logger.info(msg="write data to {}...".format(self.file_path))

        if exists(self.file_path):
            os.remove(self.file_path)

        with h5py.File(self.file_path, 'w') as f:
            f.attrs["num_points"] = self.num_points
            f.create_dataset('coord_xy', data=self.coord_xy)
            f.create_dataset('point_id', data=self.point_id)
            f.create_dataset('phase', data=self.phase)

    def open(self, path_inputs: str, other_file_path: str = None):
        """Read data from file.

        Read stored information from already existing .h5 file. This can be the file of the object itself. If the
        data should be read from another file, the path to this file can be given as 'other_file_path'. Thereby, a new
        Points object can be created with the data of another Points object.

        Parameters
        ----------
        path_inputs: str
            path to input files (slcStack.h5, geometryRadar.h5).
        other_file_path: str
            path to other .h5 file (default: None).
        """
        # 1) read own data: coord_xy, phase, point_id, num_points, reference_point_idx
        if other_file_path is not None:
            path = other_file_path
        else:
            path = self.file_path
        self.logger.info(msg="read from {}".format(path))

        with h5py.File(path, 'r') as f:
            self.num_points = f.attrs["num_points"]
            self.coord_xy = f["coord_xy"][:]
            self.point_id = f["point_id"][:]
            self.phase = f["phase"][:]

        self.openExternalData(path_inputs=path_inputs)

    def openExternalData(self, *, path_inputs: str):
        """Load data which is stored in slcStack.h5, geometryRadar.h5, ifg_network.h5 and coordinates_utm.h5."""
        # 1) read IfgNetwork
        self.ifg_net_obj.open(path=join(dirname(self.file_path), "ifg_network.h5"))

        # 2) read metadata from slcStack
        slc_stack_obj = slcStack(join(path_inputs, "slcStack.h5"))
        slc_stack_obj.open(print_msg=False)
        self.wavelength = np.float64(slc_stack_obj.metadata["WAVELENGTH"])
        self.length = slc_stack_obj.length  # y-coordinate axis (azimut)
        self.width = slc_stack_obj.width  # x-coordinate axis (range)

        # 3) read from geometry file
        mask = self.createMask()

        geom_path = join(path_inputs, "geometryRadar.h5")

        # load geometry data
        loc_inc, meta = readfile.read(geom_path, datasetName='incidenceAngle')
        loc_inc *= np.pi / 180  # in [rad]
        slant_range = readfile.read(geom_path, datasetName='slantRangeDistance')[0]
        height = readfile.read(geom_path, datasetName='height')[0]
        lat = readfile.read(geom_path, datasetName='latitude')[0]
        lon = readfile.read(geom_path, datasetName='longitude')[0]

        self.loc_inc = loc_inc[mask].ravel()
        self.slant_range = slant_range[mask].ravel()
        self.height = height[mask].ravel()
        self.coord_lalo = np.array([lat[mask].ravel(), lon[mask].ravel()]).transpose()

        # 4) read UTM coordinates
        coord_utm_obj = CoordinatesUTM(file_path=join(dirname(self.file_path), "coordinates_utm.h5"),
                                       logger=self.logger)
        coord_utm_obj.open()
        self.coord_utm = coord_utm_obj.coord_utm[:, mask].transpose()

    def createMask(self):
        """Create a mask.

        Create a mask in the size of the radar image which is used to read the geometry and SLC data for the selected
        scatterers.
        """
        mask = np.zeros((self.length, self.width), dtype=np.bool_)
        tmp = [tuple([c[0], c[1]]) for c in self.coord_xy]
        for i in tmp:
            mask[i] = True
        return mask

    def addPointsFromObj(self, *, new_point_id: np.ndarray, new_coord_xy: np.ndarray, new_phase: np.ndarray,
                         new_num_points: int, path_inputs: str):
        """Add new points and their attributes to the existing data.

        Parameters
        ----------
        new_point_id: np.ndarray
            point_id of the new scatterers.
        new_coord_xy: np.ndarray
            radar coordinates of the new scatterers.
        new_phase: np.ndarray
            phase of the new scatterers.
        new_num_points: int
            number of new points.
        path_inputs: str
            path to input files (slcStack.h5, geometryRadar.h5).
        """
        self.point_id = np.append(self.point_id, new_point_id)
        self.coord_xy = np.append(self.coord_xy, new_coord_xy, axis=0)
        self.phase = np.append(self.phase, new_phase, axis=0)
        self.num_points += new_num_points

        # all data must be ordered, so that all external data can be loaded correctly
        sort_idx = np.argsort(self.point_id)
        self.point_id = self.point_id[sort_idx]
        self.coord_xy = self.coord_xy[sort_idx, :]
        self.phase = self.phase[sort_idx, :]
        # refresh by reopening all external data
        self.openExternalData(path_inputs=path_inputs)

    def removePoints(self, mask: np.ndarray = None, *, keep_id: [np.ndarray, list], path_inputs: str):
        """Remove all entries from specified points.

        The possible options exist for removing the points:
        a) Keep all points which are set to True in a 'mask' with size (num_points x 1). Or
        b) Keep all points whose ID is listed in keep_id. The rest of the points will be removed.

        Parameters
        ----------
        mask: np.ndarray
            mask to select points to be kept, rest will be removed (default: None).
        keep_id: np.ndarray
            list of point_id to keep.
        path_inputs: str
            path to input files (slcStack.h5, geometryRadar.h5).
        """
        if mask is None:
            mask = np.ones((self.num_points,), dtype=np.bool_)
            for p in self.point_id:
                if p not in keep_id:
                    mask[self.point_id == p] = False
        self.point_id = self.point_id[mask]
        self.coord_xy = self.coord_xy[mask, :]
        self.phase = self.phase[mask, :]
        self.num_points = mask[mask].shape[0]
        # refresh by reopening all external data
        self.openExternalData(path_inputs=path_inputs)


class Network:
    """Spatial network of PS candidates."""

    def __init__(self, *, file_path: str, logger: Logger):
        """Init.

        Parameters
        ----------
        file_path: str
            absolute path to working directory for creating/loading 'psNetwork.h5'
        logger: Logger
            Logging handler.
        """
        self.num_arcs = None
        self.gamma = None
        self.arcs = None
        self.slant_range = None
        self.loc_inc = None
        self.phase = None
        self.vel = None
        self.demerr = None
        self.ifg_net_obj = None
        self.width = None
        self.length = None
        self.wavelength = None
        self.file_path = file_path
        self.logger = logger

    def writeToFile(self):
        """Write all existing data to psNetwork.h5 file."""
        self.logger.info(msg="write data to {}...".format(self.file_path))

        if exists(self.file_path):
            os.remove(self.file_path)

        with h5py.File(self.file_path, 'w') as f:
            f.attrs["num_arcs"] = self.num_arcs
            f.create_dataset('arcs', data=self.arcs)
            f.create_dataset('phase', data=self.phase)
            f.create_dataset('loc_inc', data=self.loc_inc)
            f.create_dataset('slant_range', data=self.slant_range)

    def open(self, *, path_inputs: str):
        """Read stored information from existing .h5 file."""
        with h5py.File(self.file_path, 'r') as f:
            self.num_arcs = f.attrs["num_arcs"]
            self.arcs = f["arcs"][:]
            self.phase = f["phase"][:]
            self.loc_inc = f["loc_inc"][:]
            self.slant_range = f["slant_range"][:]
        self.openExternalData(path_inputs=path_inputs)

    def openExternalData(self, *, path_inputs: str):
        """Read data from slcStack.h5 and IfgNetwork.h5 files."""
        slc_stack_obj = slcStack(join(path_inputs, "slcStack.h5"))
        slc_stack_obj.open(print_msg=False)
        self.wavelength = np.float64(slc_stack_obj.metadata["WAVELENGTH"])
        self.length = slc_stack_obj.length  # y-coordinate axis (azimut)
        self.width = slc_stack_obj.width  # x-coordinate axis (range)

        # 3) read IfgNetwork
        self.ifg_net_obj = IfgNetwork()
        self.ifg_net_obj.open(path=join(dirname(self.file_path), "ifg_network.h5"))

    def computeArcObservations(self, *, point_obj: Points, arcs: np.ndarray):
        """Compute the phase observations for each arc.

        Compute double difference phase observations, i.e. the phase differences for each arc in the network from the
        phase of the two scatterers connected by the arc.

        Parameters
        ----------
        point_obj: Points
            object of class Points.
        arcs: np.ndarray
            Array with the indices of the points connected by an arc.
        """
        self.arcs = arcs
        self.num_arcs = self.arcs.shape[0]
        self.logger.info(msg="no. arcs:\t{}".format(self.num_arcs))

        self.phase = np.zeros((self.num_arcs, point_obj.ifg_net_obj.num_ifgs))
        self.loc_inc = np.zeros((self.num_arcs,))
        self.slant_range = np.zeros((self.num_arcs,))
        for idx, arc in enumerate(self.arcs):
            self.phase[idx, :] = np.angle(
                np.exp(1j * point_obj.phase[arc[0], :]) * np.conjugate(np.exp(1j * point_obj.phase[arc[1], :])))
            self.loc_inc[idx] = np.mean([point_obj.loc_inc[arc[0]], point_obj.loc_inc[arc[1]]])
            self.slant_range[idx] = np.mean([point_obj.slant_range[arc[0]], point_obj.slant_range[arc[1]]])

        self.logger.info(msg="ifg arc observations created.")

    def removeArcs(self, *, mask: np.ndarray):
        """Remove arcs from the list of arcs in the network.

        Parameter
        ---------
        mask: np.ndarray
            mask to select arcs to be kept, rest will be removed.
        """
        self.demerr = self.demerr[mask]
        self.vel = self.vel[mask]
        self.phase = self.phase[mask, :]
        self.loc_inc = self.loc_inc[mask]
        self.slant_range = self.slant_range[mask]
        self.arcs = np.array(self.arcs)
        self.arcs = self.arcs[mask, :]
        self.gamma = self.gamma[mask]
        self.num_arcs = mask[mask].shape[0]


class NetworkParameter(Network):
    """Spatial Network with the estimated parameters of each arc in the network."""

    def __init__(self, *, file_path: str, logger: Logger):
        """Init."""
        super().__init__(file_path=file_path, logger=logger)
        self.gamma = None
        self.vel = None
        self.demerr = None
        self.slant_range = None
        self.loc_inc = None
        self.phase = None
        self.arcs = None
        self.num_arcs = None
        self.logger = logger

    def prepare(self, *, net_obj: Network, demerr: np.ndarray, vel: np.ndarray, gamma: np.ndarray):
        """Prepare.

        Parameter
        -----------
        net_obj: Network
            object of class Network.
        demerr: np.ndarray
            estimated DEM error for each arc in the network.
        vel: np.ndarray
            estimated velocity for each arc in the network.
        gamma: np.ndarray
            estimated temporal coherence for each arc in the network.
        """
        self.num_arcs = net_obj.num_arcs
        self.arcs = net_obj.arcs
        self.phase = net_obj.phase
        self.loc_inc = net_obj.loc_inc
        self.slant_range = net_obj.slant_range
        self.demerr = demerr
        self.vel = vel
        self.gamma = gamma

    def writeToFile(self):
        """Write DEM error, velocity and temporal coherence to file."""
        super().writeToFile()

        with h5py.File(self.file_path, 'r+') as f:  # append existing file
            f.create_dataset('demerr', data=self.demerr)
            f.create_dataset('vel', data=self.vel)
            f.create_dataset('gamma', data=self.gamma)

    def open(self, *, path_inputs: str):
        """Read data from file."""
        super().open(path_inputs=path_inputs)

        with h5py.File(self.file_path, 'r') as f:
            self.demerr = f["demerr"][:]
            self.vel = f["vel"][:]
            self.gamma = f["gamma"][:]
