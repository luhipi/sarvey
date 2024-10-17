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

"""Module for correcting the geolocation of the scatterers."""
import logging
from os.path import join
import numpy as np

from miaplpy.objects.slcStack import slcStack

from sarvey.objects import Points


def getHeading(input_path: str, logger: logging.Logger):
    """
    Read heading angle from slcStack.h5.

    Parameters
    ----------
    input_path: str
        Path to directory containing 'slcStack.h5' and 'geometryRadar.h5'.
    logger: Logger
        Logger handle

    Returns
    -------
    heading_angle: float
        heading angle of the satellite in radians
        for ascending ~ -12*pi/180
        for descending ~ 190*pi/180
    """
    # get heading from slcStack.h5
    slc_stack_file = join(input_path, 'slcStack.h5')
    slc_stack_obj = slcStack(slc_stack_file)
    try:
        meta_dict = slc_stack_obj.get_metadata()
        lower_case_meta_dict = {k.lower(): v for k, v in meta_dict.items()}

        heading_angle = float(lower_case_meta_dict["heading"])
        logger.info(msg=f"Heading_angle of satellite: {heading_angle} deg")
        heading_angle = np.deg2rad(heading_angle)
    except Exception as exc:
        logger.error(f'Failed to retrieve heading angle from {slc_stack_file}: {exc}')
        raise Exception
    return heading_angle


def calculateGeolocationCorrection(*, path_geom: str, point_obj: Points, demerr: np.array, logger: logging.Logger):
    """
    Calculate geolocation correction.

    Parameters
    ----------
    path_geom: str
        Path to directory containing 'slcStack.h5' or 'geometryRadar.h5'.
    point_obj: Points
        Point object with incidence angle for points
    demerr: np.array
        Array of dem error per pixel
    logger: Logger
        Logger handle

    Returns
    -------
    coord_correction: np.array
        array of geolocation corrections, two columns [x_correction, y_correction] per point.
    """
    heading_angle = getHeading(input_path=path_geom, logger=logger)

    coord_correction = np.zeros_like(point_obj.coord_xy, dtype=float)
    coord_correction[:, 0] = demerr * np.cos(heading_angle) / np.tan(point_obj.loc_inc)
    coord_correction[:, 1] = -demerr * np.sin(heading_angle) / np.tan(point_obj.loc_inc)

    return coord_correction
