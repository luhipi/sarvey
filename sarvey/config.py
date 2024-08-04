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

"""Configuration module for SARvey."""
import os
import json
from datetime import date
from json import JSONDecodeError
from typing import Optional
from pydantic import BaseModel, Field, validator, Extra


class DataDirectories(BaseModel, extra=Extra.forbid):
    """Template for settings in config file."""

    path_inputs: str = Field(
        title="The path to the input data directory.",
        description="Set the path of the input data directory.",
        default="inputs/"
    )

    path_outputs: str = Field(
        title="The path to the processing output data directory.",
        description="Set the path of the processing output data directory.",
        default="outputs/"
    )

    @validator('path_inputs')
    def checkPathInputs(cls, v):
        """Check if the input path exists."""
        if v == "":
            raise ValueError("Empty string is not allowed.")
        if not os.path.exists(os.path.abspath(v)):
            raise ValueError(f"path_inputs is invalid: {os.path.abspath(v)}")
        if not os.path.exists(os.path.join(os.path.abspath(v), "slcStack.h5")):
            raise ValueError(f"'slcStack.h5' does not exist: {v}")
        if not os.path.exists(os.path.join(os.path.abspath(v), "geometryRadar.h5")):
            raise ValueError(f"'geometryRadar.h5' does not exist: {v}")
        return v


class Processing(BaseModel, extra=Extra.forbid):
    """Template for settings in config file."""

    num_cores: int = Field(
        title="Number of cores",
        description="Set the number of cores for parallel processing.",
        default=50
    )

    num_patches: int = Field(
        title="Number of patches",
        description="Set the number of patches for processing large areas patch-wise.",
        default=1
    )

    temporal_unwrapping: bool = Field(
        title="Apply temporal unwrapping",
        description="Apply temporal unwrapping additionally to spatial unwrapping.",
        default=True
    )

    unwrapping_method: str = Field(
        title="Spatial unwrapping method",
        description="Select spatial unwrapping method from 'ilp' and 'puma'.",
        default='puma'
    )

    @validator('num_cores')
    def checkNumCores(cls, v):
        """Check if the number of cores is valid."""
        if v <= 0:
            raise ValueError("Number of cores must be greater than zero.")
        return v

    @validator('num_patches')
    def checkNumPatches(cls, v):
        """Check if the number of patches is valid."""
        if v <= 0:
            raise ValueError("Number of patches must be greater than zero.")
        return v

    @validator('unwrapping_method')
    def checkUnwMethod(cls, v):
        """Check if unwrapping_method is valid."""
        if (v != "ilp") & (v != "puma"):
            raise ValueError("Unwrapping method must be either 'ilp' or 'puma'.")
        return v


class PhaseLinking(BaseModel, extra=Extra.forbid):
    """Template for settings in config file."""

    phase_linking: bool = Field(
        title="Use phase linking results.",
        description="Use pixels selected in phase linking.",
        default=False
    )

    path_inverted: str = Field(
        title="The path to the phase linking inverted data directory.",
        description="Set the path of the inverted data directory.",
        default="inverted/"
    )

    num_siblings: int = Field(
        title="Sibling threshold.",
        description="Threshold on the number of siblings applied during phase linking to distinguish PS from DS"
                    "candidates.",
        default=20
    )

    spatial_mask_file_pl: Optional[str] = Field(
        title="Path to spatial mask file for phase linking results.",
        description="Path to the mask file, e.g. created by sarvey_mask.",
        default=""
    )

    use_ps: bool = Field(
        title="Use point-like scatterers.",
        description="Use point-like scatterers (pixels with a low number of siblings) selected in phase linking."
                    "Is applied, only if 'phase_linking' is true.",
        default=False
    )

    path_mask_file_ps: str = Field(
        title="The path to the mask file for ps pixels from phase linking.",
        description="Set the path of the 'maskPS.h5' file (optional).",
        default="maskPS.h5"
    )

    @validator('path_inverted')
    def checkPathInverted(cls, v, values):
        """Check if the inverted path exists."""
        if values["phase_linking"]:
            if v == "":
                raise ValueError("Empty string is not allowed.")
            if not os.path.exists(os.path.abspath(v)):
                raise ValueError(f"path_inverted is invalid: {os.path.abspath(v)}")
            if not os.path.exists(os.path.join(os.path.abspath(v), "phase_series.h5")):
                raise ValueError(f"'phase_series.h5' does not exist: {v}")
        return v

    @validator('num_siblings')
    def checkNumSiblings(cls, v, values):
        """Check is no_siblings is valid."""
        if not values["phase_linking"]:
            if v < 1:
                raise ValueError("'num_siblings' has to be greater than 0.")
        return v

    @validator('spatial_mask_file_pl')
    def checkSpatialMaskPath(cls, v, values):
        """Check if the path is correct."""
        if values["phase_linking"]:
            if v == "" or v is None:
                return None
            else:
                if not os.path.exists(os.path.abspath(v)):
                    raise ValueError(f"spatial_mask_file_pl path is invalid: {v}")
        return v

    @validator('use_ps')
    def checkUsePS(cls, v, values):
        """Check if use_ps will be applied."""
        if (not values["phase_linking"]) and v:
            raise ValueError("'use_ps' will not be applied, because 'phase_linking' is set to False.")
        return v

    @validator('path_mask_file_ps')
    def checkPathMaskFilePS(cls, v, values):
        """Check if the mask file exists."""
        if values["phase_linking"] and values["use_ps"]:
            if v == "":
                raise ValueError("Empty string is not allowed.")
            if not os.path.exists(os.path.abspath(v)):
                raise ValueError(f"path_mask_file_ps is invalid: {os.path.abspath(v)}")
        return v


class Preparation(BaseModel, extra=Extra.forbid):
    """Template for settings in config file."""

    start_date: Optional[str] = Field(
        title="Start date",
        description="Format: YYYY-MM-DD.",
        default=None
    )

    stop_date: Optional[str] = Field(
        title="Stop date",
        description="Format: YYYY-MM-DD.",
        default=None
    )

    network_type: str = Field(
        title="Interferogram network type.",
        description="Set the intererogram network type: 'sb' (small baseline), 'stb' (small temporal baseline), "
                    "'stb_year' (small temporal baseline and yearly ifgs), 'delaunay' (delaunay network), "
                    "or 'star' (single-reference).",
        default="sb"
    )

    num_ifgs: Optional[int] = Field(
        title="Number of interferograms",
        description="Set the number of interferograms per image. Might be violated .",
        default=3
    )

    max_tbase: Optional[int] = Field(
        title="Maximum temporal baseline [days]",
        description="Set the maximum temporal baseline for the ifg network. (required for: 'sb')",
        default=100
    )

    filter_wdw_size: int = Field(
        title="Size of filtering window [pixel]",
        description="Set the size of window for lowpass filtering.",
        default=9
    )

    @validator('start_date', 'stop_date')
    def checkDates(cls, v):
        """Check if date format is valid."""
        if v == "":
            v = None

        if v is not None:
            try:
                date.fromisoformat(v)
            except Exception as e:
                raise ValueError(f"Date needs to be in format: YYYY-MM-DD. {e}")
        return v

    @validator('network_type')
    def checkNetworkType(cls, v):
        """Check if the ifg network type is valid."""
        if (v != "sb") and (v != "star") and (v != "delaunay") and (v != "stb") and (v != "stb_year"):
            raise ValueError("Interferogram network type has to be 'sb', 'stb', Ästb_year', 'delaunay' or 'star'.")
        return v

    @validator('num_ifgs')
    def checkNumIfgs(cls, v):
        """Check if the number of ifgs is valid."""
        if v is not None:
            if v <= 0:
                raise ValueError("Number of ifgs must be greater than zero.")
        return v

    @validator('max_tbase')
    def checkMaxTBase(cls, v):
        """Check if the value for maximum time baseline is valid."""
        if v is not None:
            if v <= 0:
                raise ValueError("Maximum baseline must be greater than zero.")
        return v

    @validator('filter_wdw_size')
    def checkFilterWdwSize(cls, v):
        """Check if the filter window size is valid."""
        if v <= 0:
            raise ValueError("Filter window size must be greater than zero.")
        return v


class ConsistencyCheck(BaseModel, extra=Extra.forbid):
    """Template for settings in config file."""

    coherence_p1: float = Field(
        title="Temporal coherence threshold for first-order points",
        description="Set the temporal coherence threshold of first-order points for the consistency check.",
        default=0.9
    )

    grid_size: int = Field(
        title="Grid size [m]",
        description="Set the grid size in [m] for the consistency check. No grid is applied if 'grid_size' is Zero.",
        default=200
    )

    spatial_mask_file_p1: Optional[str] = Field(
        title="Path to mask file",
        description="Set the path to the mask file in .h5 format.",
        default=""
    )

    knn: int = Field(
        title="Number of nearest neighbours",
        description="Set number of nearest neighbours for creating arcs.",
        default=30
    )

    max_arc_length: Optional[int] = Field(
        title="Maximum length of arcs [m]",
        description="Set the maximum length of arcs.",
        default=None
    )

    velocity_bound: float = Field(
        title="Bounds on mean velocity for temporal unwrapping [m/year]",
        description="Set the bound (symmetric) for the mean velocity estimation in temporal unwrapping.",
        default=0.1
    )

    dem_error_bound: float = Field(
        title="Bounds on DEM error for temporal unwrapping [m]",
        description="Set the bound (symmetric) for the DEM error estimation in temporal unwrapping.",
        default=100.0
    )

    num_samples: int = Field(
        title="Number of samples in the search space for temporal unwrapping",
        description="Set the number of samples evaluated along the search space for temporal unwrapping.",
        default=100
    )

    arc_coherence: float = Field(
        title="Arc coherence threshold",
        description="Set the arc coherence threshold for the consistency check.",
        default=0.6
    )

    min_num_arc: int = Field(
        title="Minimum number of arcs per point",
        description="Set the minimum number of arcs per point.",
        default=3
    )

    @validator('coherence_p1')
    def checkCoherenceP1(cls, v):
        """Check if the temporal coherence threshold is valid."""
        if v < 0:
            raise ValueError("Temporal coherence threshold cannot be negative.")
        if v > 1:
            raise ValueError("Temporal coherence threshold cannot be greater than 1.")
        return v

    @validator('grid_size')
    def checkGridSize(cls, v):
        """Check if the grid size is valid."""
        if v < 0:
            raise ValueError('Grid size cannot be negative.')
        if v == 0:
            v = None
        return v

    @validator('spatial_mask_file_p1')
    def checkSpatialMaskPath(cls, v):
        """Check if the path is correct."""
        if v == "" or v is None:
            return None
        else:
            if not os.path.exists(os.path.abspath(v)):
                raise ValueError(f"spatial_mask_file_p1 path is invalid: {v}")
            return v

    @validator('knn')
    def checkKNN(cls, v):
        """Check if the k-nearest neighbours is valid."""
        if v <= 0:
            raise ValueError('K-nearest neighbours cannot be negative or zero.')
        return v

    @validator('max_arc_length')
    def checkMaxArcLength(cls, v):
        """Check if the maximum length of arcs is valid."""
        if v is None:
            return 999999
        if v <= 0:
            raise ValueError('Maximum arc length must be positive.')
        return v

    @validator('velocity_bound')
    def checkVelocityBound(cls, v):
        """Check if the velocity bound is valid."""
        if v <= 0:
            raise ValueError('Velocity bound cannot be negative or zero.')
        return v

    @validator('dem_error_bound')
    def checkDEMErrorBound(cls, v):
        """Check if the DEM error bound is valid."""
        if v <= 0:
            raise ValueError('DEM error bound cannot be negative or zero.')
        return v

    @validator('num_samples')
    def checkNumSamples(cls, v):
        """Check if the number of samples for the search space is valid."""
        if v <= 0:
            raise ValueError('Number of samples cannot be negative or zero.')
        return v

    @validator('arc_coherence')
    def checkArcCoherence(cls, v):
        """Check if the arc coherence threshold is valid."""
        if v < 0:
            raise ValueError('Arc coherence threshold cannot be negativ.')
        if v > 1:
            raise ValueError('Arc coherence threshold cannot be greater than 1.')
        return v

    @validator('min_num_arc')
    def checkMinNumArc(cls, v):
        """Check if the minimum number of arcs is valid."""
        if v < 0:
            raise ValueError('Velocity bound cannot be negative.')
        return v


class Unwrapping(BaseModel, extra=Extra.forbid):
    """Template for settings in config file."""

    knn: int = Field(
        title="Number of nearest neighbours",
        description="Set number of nearest neighbours for spatial unwrapping with PUMA. Not used, if"
                    "'use_temporal_unwrapping_arcs' is 'true'.",
        default=1
    )

    use_temporal_unwrapping_arcs: bool = Field(
        title="Use arcs from temporal unwrapping",
        description="If true, use same arcs from temporal unwrapping. If false, apply new delaunay and knn"
                    "triangulation.",
        default=True
    )

    @validator('knn')
    def checkKNN(cls, v):
        """Check if the k-nearest neighbours is valid."""
        if v <= 0:
            raise ValueError('K-nearest neighbours cannot be negative or zero.')
        return v


class Filtering(BaseModel, extra=Extra.forbid):
    """Template for filtering settings in config file."""

    skip_filtering: bool = Field(
        title="Skip filtering step.",
        description="Set whether to skip filtering step.",
        default=False
    )

    interpolation_method: str = Field(
        title="Spatial interpolation method.",
        description="Method for interpolating atmosphere in space ('linear', 'cubic' or 'kriging').",
        default="kriging"
    )

    grid_size: int = Field(
        title="Grid size [m].",
        description="Set the grid size for spatial filtering.",
        default=1000
    )

    spatial_mask_file_p2: Optional[str] = Field(
        title="Path to spatial mask file.",
        description="Path to the mask file, e.g. created by sarvey_mask.",
        default=""
    )

    use_moving_points: bool = Field(
        title="Use moving points",
        description="Set whether to use moving points in the filtering step.",
        default=True
    )

    max_auto_corr: float = Field(
        title="Max auto correlation.",
        description="Set temporal autocorrelation threshold for the selection of stable/linearly moving points.",
        default=0.3
    )

    @validator('interpolation_method')
    def checkInterpolationMethod(cls, v):
        """Check if the interpolation method is valid."""
        if (v.lower() != "linear") and (v.lower() != "cubic") and (v.lower() != "kriging"):
            raise ValueError("Method for interpolating atmosphere in space needs to be either 'linear', 'cubic' "
                             "or 'kriging'.")
        return v

    @validator('grid_size')
    def checkGridSize(cls, v):
        """Check if the grid size is valid."""
        if v < 0:
            raise ValueError("Grid size cannot be negative.")
        else:
            return v

    @validator('spatial_mask_file_p2')
    def checkSpatialMaskPath(cls, v):
        """Check if the path is correct."""
        if v == "" or v is None:
            return None
        else:
            if not os.path.exists(os.path.abspath(v)):
                raise ValueError(f"spatial_mask_file_p2 path is invalid: {v}")
        return v

    @validator('max_auto_corr')
    def checkMaxAutoCorr(cls, v):
        """Check if the value is correct."""
        if v < 0 or v > 1:
            raise ValueError(f"max_auto_corr is not between 0 and 1: {v}")
        return v


class Densification(BaseModel, extra=Extra.forbid):
    """Template for densification settings in config file."""

    coherence_p2: float = Field(
        title="Temporal coherence threshold",
        description="Set the temporal coherence threshold for the densification step.",
        default=0.8
    )

    coherence_threshold: float = Field(
        title="Coherence threshold for densification",
        description="Set coherence threshold for densification.",
        default=0.5
    )

    num_connections_p1: int = Field(
        title="Number of connections in temporal unwrapping.",
        description="Set number of connections between second-order point and closest first-order points for temporal "
                    "unwrapping.",
        default=5
    )

    num_connections_p2: int = Field(
        title="Number of connections in consistency check with neighbouring points.",
        description="Set number of connections between unwrapped second-order point and closest second-order points for"
                    " temporal checking consistency.",
        default=10
    )

    max_distance_p1: int = Field(
        title="Maximum distance to nearest first-order point [m]",
        description="Set threshold on the distance between first-order points and to be temporally unwrapped"
                    "second-order point.",
        default=2000
    )

    velocity_bound: float = Field(
        title="Bounds on mean velocity for temporal unwrapping [m/year]",
        description="Set the bound (symmetric) for the mean velocity in temporal unwrapping.",
        default=0.15
    )

    dem_error_bound: float = Field(
        title="Bounds on DEM error for temporal unwrapping [m]",
        description="Set the bound (symmetric) for the DEM error estimation in temporal unwrapping.",
        default=100.0
    )

    num_samples: int = Field(
        title="Number of samples in the search space for temporal unwrapping",
        description="Set the number of samples evaluated along the search space for temporal unwrapping.",
        default=100
    )

    knn: int = Field(
        title="Number of nearest neighbours",
        description="Set number of nearest neighbours for creating arcs.",
        default=1
    )

    @validator('coherence_p2')
    def checkTempCohThrsh2(cls, v):
        """Check if the temporal coherence threshold is valid."""
        if v < 0:
            raise ValueError("Temporal coherence threshold cannot be negative.")
        if v > 1:
            raise ValueError("Temporal coherence threshold cannot be greater than 1.")
        return v

    @validator('coherence_threshold')
    def checkCoherenceThresh(cls, v):
        """Check if coherence_threshold is valid."""
        if v < 0 or v > 1:
            raise ValueError(f"coherence_threshold is not between 0 and 1: {v}")
        return v

    @validator('num_connections_p1')
    def checkNumConn1(cls, v):
        """Check if num_connections_p1 are valid."""
        if v <= 0:
            raise ValueError(f"num_connections_p1 must be greater than 0: {v}")
        return v

    @validator('num_connections_p1')
    def checkNumConn2(cls, v):
        """Check if num_connections_p2 are valid."""
        if v < 0:
            raise ValueError(f"num_connections_p2 cannot be negative: {v}")
        return v

    @validator('max_distance_p1')
    def checkMaxDistanceP1(cls, v):
        """Check if the maximum distance to nearest first-order points is valid."""
        if v < 0:
            raise ValueError('Maximum distance to first-order points cannot be negative.')
        return v

    @validator('velocity_bound')
    def checkVelocityBound(cls, v):
        """Check if the velocity bound is valid."""
        if v <= 0:
            raise ValueError('Velocity bound cannot be negative or zero.')
        return v

    @validator('dem_error_bound')
    def checkDEMErrorBound(cls, v):
        """Check if the DEM error bound is valid."""
        if v <= 0:
            raise ValueError('DEM error bound cannot be negative or zero.')
        return v

    @validator('num_samples')
    def checkNumSamples(cls, v):
        """Check if the number of samples for the search space is valid."""
        if v <= 0:
            raise ValueError('Number of samples cannot be negative or zero.')
        return v

    @validator('knn')
    def checkKNN(cls, v):
        """Check if the k-nearest neighbours is valid."""
        if v <= 0:
            raise ValueError('K-nearest neighbours cannot be negative or zero.')
        return v


class Logging(BaseModel):
    """Template for logger settings."""

    logging_level: str = Field(
        title="Logging level.",
        description="Set loggig level.",
        default="INFO"
    )

    logfile_path: str = Field(
        title="Logfile Path.",
        description="Path to directory where the logfiles should be saved.",
        default="logfiles/"
    )

    @validator('logging_level')
    def checkLoggingLevel(cls, v):
        """Check if the logging level is valid."""
        if v == "":
            raise ValueError("Empty string is not allowed.")
        v = v.upper()
        if v not in ["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]:
            raise ValueError("Logging level must be one of ('CRITICAL', 'ERROR', "
                             "'WARNING', 'INFO', 'DEBUG', 'NOTSET').")
        return v


class Config(BaseModel):
    """Configuration for SAR4Infra sarvey."""

    # title has to be the name of the class. Needed for creating default file
    data_directories: DataDirectories = Field(
        title="DataDirectories", description=""
    )

    logging: Logging = Field(
        title="Logging", description=""
    )

    processing: Processing = Field(
        title="Processing", description=""
    )

    phase_linking: PhaseLinking = Field(
        title="PhaseLinking", description=""
    )

    preparation: Preparation = Field(
        title="Preparation", description=""
    )

    consistency_check: ConsistencyCheck = Field(
        title="ConsistencyCheck", description=""
    )

    unwrapping: Unwrapping = Field(
        title="Unwrapping", description=""
    )

    filtering: Filtering = Field(
        title="Filtering", description=""
    )

    densification: Densification = Field(
        title="Densification", description=""
    )


def loadConfiguration(*, path: str) -> dict:
    """Load configuration json file.

    Parameters
    ----------
    path : str
        Path to the configuration json file.

    Returns
    -------
    : dict
        A dictionary containing configurations.

    Raises
    ------
    JSONDecodeError
        If failed to parse the json file to the dictionary.
    FileNotFoundError
        Config file not found.
    IOError
        Invalid JSON file.
    ValueError
        Invalid value for configuration object.
    """
    try:
        with open(path) as config_fp:
            config = json.load(config_fp)
            config = Config(**config).dict(by_alias=True)
    except JSONDecodeError as e:
        raise IOError(f'Failed to load the configuration json file => {e}')
    return config
