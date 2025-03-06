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

"""Utils module for SARvey."""
import time
from os.path import exists, join

import numpy as np
from scipy.sparse.linalg import lsqr
from typing import Union
from logging import Logger

from mintpy.utils import ptime

from sarvey.objects import Points, NetworkParameter, Network, BaseStack, AmplitudeImage
from sarvey.ifg_network import IfgNetwork


def convertBboxToBlock(*, bbox: tuple):
    """ConvertBboxToBlock. read box and write2hdf5_block block have different definitions."""
    block = None
    if len(bbox) == 4:
        block = (bbox[1], bbox[3], bbox[0], bbox[2])
    if len(bbox) == 6:
        block = (bbox[2], bbox[5], bbox[1], bbox[4], bbox[0], bbox[3])
    return block


def invertIfgNetwork(*, phase: np.ndarray, num_points: int, ifg_net_obj: IfgNetwork, logger: Logger,
                     lifetime_ifgs: np.ndarray = None, lifetime_images: np.ndarray = None):
    """Wrap the ifg network inversion running in parallel.

    This function does not use multiprocessing. Matrix inversion in scipy is already parallelized.

    Parameters
    ----------
    phase: np.ndarray
        interferometric phases of the points.
    num_points: int
        number of points.
    ifg_net_obj: IfgNetwork
        instance of class IfgNetwork.
    logger: Logger
        logging handler
    lifetime_ifgs: np.ndarray
        lifetime of each point for each interferogram (optional).
    lifetime_images: np.ndarray
        lifetime of each point for each image (optional).

    Returns
    -------
        phase_ts: np.ndarray
            inverted phase time series of the points.
    """
    msg = "#" * 10
    msg += " INVERT IFG NETWORK "
    msg += "#" * 10
    logger.info(msg=msg)

    # for legacy, if only continuously coherent scatterers are processed:
    if lifetime_ifgs is None:
        lifetime_ifgs = np.ones((num_points, ifg_net_obj.num_ifgs), dtype=np.bool_)
    if lifetime_images is None:
        lifetime_images = np.ones((num_points, ifg_net_obj.num_images), dtype=np.bool_)

    start_time = time.time()
    design_mat = ifg_net_obj.getDesignMatrix()

    phase_ts = np.zeros((num_points, ifg_net_obj.num_images), dtype=np.float32)
    phase_ts[~lifetime_images] = np.nan

    prog_bar = ptime.progressBar(maxValue=num_points)
    for p in range(num_points):
        mask_coherent = lifetime_images[p, :].copy()
        ref_idx = np.where(mask_coherent)[0][0]  # first coherent images of the point
        mask_coherent[ref_idx] = False

        design_mat_adjusted = design_mat[lifetime_ifgs[p, :], :]
        design_mat_adjusted = design_mat_adjusted[:, lifetime_images[p, :]]
        # remove reference date, caution: reference date is always the first date inside the cropped design matrix.
        design_mat_adjusted = np.delete(design_mat_adjusted, 0, axis=1)

        phase_ts[p, mask_coherent] = lsqr(design_mat_adjusted, phase[p, lifetime_ifgs[p, :]])[0]
        prog_bar.update(value=p + 1, every=np.ceil(num_points / 100), suffix='{}/{} points'.format(p + 1, num_points))

    m, s = divmod(time.time() - start_time, 60)
    logger.debug(msg='time used: {:02.0f} mins {:02.1f} secs.'.format(m, s))
    return phase_ts


def predictPhase(*, obj: [NetworkParameter, Points], vel: np.ndarray = None, demerr: np.ndarray = None,
                 ifg_space: bool = True, logger: Logger):
    """Predicts the phase time series based on the estimated parameters DEM error and mean velocity.

    Can be used for both arc phase or point phase. Wrapper function for 'predictPhaseCore(...)'

    Parameters
    ----------
    obj: Union[NetworkParameter, Points]
        object of either 'networkParameter' or 'points'. If instance of 'points' is given, 'vel' and 'demerr'
        also need to be specified.
    vel: np.ndarray
        velocity for each sample (default: None)
    demerr: np.ndarray
        dem error for each sample (default: None).
    ifg_space: bool
        set to True if the phase shall be predicted in interferogram space. If False, phase will be predicted
        in acquisition space. (default: True)
    logger: Logger
        Logging handler.

    Returns
    -------
        pred_phase_demerr: np.ndarray
            predicted phase from DEM error
        pred_phase_vel: np.ndarray
            predicted phase from velocity

    Raises
    ------
    ValueError
        vel or demerr is none
    TypeError
        obj is of the wrong type
    """
    if isinstance(obj, Points):
        if (vel is None) or (demerr is None):
            logger.error(msg="Both 'vel' and 'demerr' are needed if 'obj' is instance of class 'points'!")
            raise ValueError
        pred_phase_demerr, pred_phase_vel = predictPhaseCore(
            ifg_net_obj=obj.ifg_net_obj,
            wavelength=obj.wavelength,
            vel=vel,
            demerr=demerr,
            slant_range=obj.slant_range,
            loc_inc=obj.loc_inc,
            ifg_space=ifg_space
        )
    elif isinstance(obj, NetworkParameter):
        pred_phase_demerr, pred_phase_vel = predictPhaseCore(
            ifg_net_obj=obj.ifg_net_obj,
            wavelength=obj.wavelength,
            vel=obj.vel,
            demerr=obj.demerr,
            slant_range=obj.slant_range,
            loc_inc=obj.loc_inc,
            ifg_space=ifg_space
        )
    else:
        logger.error(msg="'obj' must be instance of 'points' or 'networkParameter'!")
        raise TypeError
    return pred_phase_demerr, pred_phase_vel


def predictPhaseCore(*, ifg_net_obj: IfgNetwork, wavelength: float, vel: np.ndarray,
                     demerr: np.ndarray, slant_range: np.ndarray, loc_inc: np.ndarray, ifg_space: bool = True):
    """Predicts the phase time series based on the estimated parameters DEM error and mean velocity.

    Can be used for both arc phase or point phase.

    Parameters
    ----------
    ifg_net_obj: IfgNetwork
        instance of class ifgNetwork
    wavelength: float
        wavelength in [m]
    vel: np.ndarray
        velocity for each sample
    demerr: np.ndarray
        dem error for each sample
    slant_range: np.ndarray
        slant range distance for each sample
    loc_inc: np.ndarray
        local incidence angle for each sample
    ifg_space: bool
        set to True if the phase shall be predicted in interferogram space. If False, phase will be
        predicted in acquisition space. (default: True)

    Returns
    -------
        pred_phase_demerr: np.ndarray
            predicted phase from DEM error
        pred_phase_vel: np.ndarray
            predicted phase from velocity
    """
    factor = 4 * np.pi / wavelength

    if ifg_space:
        tbase = ifg_net_obj.tbase_ifg
        pbase = ifg_net_obj.pbase_ifg
    else:
        tbase = ifg_net_obj.tbase
        pbase = ifg_net_obj.pbase

    # compute phase due to DEM error
    pred_phase_demerr = factor * pbase[:, np.newaxis] / (slant_range * np.sin(loc_inc))[np.newaxis, :] * demerr

    # compute phase due to velocity
    pred_phase_vel = factor * tbase[:, np.newaxis] * vel

    return pred_phase_demerr.T, pred_phase_vel.T


def predictPhaseSingle(*, demerr: float, vel: float, slant_range: float, loc_inc: float,
                       ifg_net_obj: IfgNetwork, wavelength: float, only_vel: bool = False, ifg_space: bool = True):
    """Predict the phase time series for only one point based on the estimated parameters DEM error and mean velocity.

    Can be used for both arc phase or point phase.

    Parameters
    ----------
    demerr: float
        DEM error (scalar)
    vel: float
        mean velocity (scalar)
    slant_range: float
        slant range distance in [m] (scalar)
    loc_inc: float
        local incidence angle in [rad] (scalar)
    ifg_net_obj: IfgNetwork
        object of class IfgNetwork
    wavelength: float
        radar wavelength in [m]
    only_vel: bool
        set to True if only the mean velocity shall be predicted (default: False)
    ifg_space: bool
        set to True if the phase shall be predicted in interferogram space. If False, phase will be predicted in
        acquisition space. (default: True)

    Returns
    -------
        pred_phase: np.ndarray
            predicted phase
    """
    factor = 4 * np.pi / wavelength

    if ifg_space:
        tbase = ifg_net_obj.tbase_ifg
        pbase = ifg_net_obj.pbase_ifg
        num_time = ifg_net_obj.num_ifgs
    else:
        tbase = ifg_net_obj.tbase
        pbase = ifg_net_obj.pbase
        num_time = ifg_net_obj.num_images

    if only_vel:
        a = np.zeros((num_time, 1))
    else:
        a = np.zeros((num_time, 2))
    a[:, 0] = factor * tbase

    if only_vel:
        pred_phase = np.matmul(a, np.array([vel])).reshape((-1,))
    else:
        a[:, 1] = factor * pbase / (slant_range * np.sin(loc_inc))
        pred_phase = np.matmul(a, np.array([vel, demerr])).reshape((-1,))

    return pred_phase


def estimateParameters(*, obj: Union[Points, Network], estimate_ref_atmo: bool = True, ifg_space: bool = True):
    """Estimate the parameters either per point or per arc.

    Parameters are velocity and DEM error (or additionally reference APS).

    Parameters
    ----------
    obj: Union[Points, Network]
        object of either network, networkParameter, points or pointsParameters
    estimate_ref_atmo: bool
        set to True if APS of reference date shall be estimated. corresponds to offset of linear
        motion model (default: False).
    ifg_space: bool
        set to True if the phase shall be predicted in interferogram space. If False, phase will be
        predicted in acquisition space. (default: True)

    Returns
    -------
    vel: np.ndarray
        velocity for each point
    demerr: np.ndarray
        dem error for each point
    ref_atmo: np.ndarray
        reference APS for each point
    omega:
        sum of squared residuals
    v_hat:
        residuals
    """
    num = obj.phase.shape[0]  # either number of points or number of arcs

    if ifg_space:
        tbase = obj.ifg_net_obj.tbase_ifg
        pbase = obj.ifg_net_obj.pbase_ifg
        num_time = obj.ifg_net_obj.num_ifgs
    else:
        tbase = obj.ifg_net_obj.tbase
        pbase = obj.ifg_net_obj.pbase
        num_time = obj.ifg_net_obj.num_images

    vel = np.zeros((num,), dtype=np.float32)
    demerr = np.zeros((num,), dtype=np.float32)
    omega = np.zeros((num,), dtype=np.float32)
    coherence = np.zeros((num,), dtype=np.float32)
    v_hat = np.zeros((num, num_time), dtype=np.float32)

    ref_atmo = None
    if estimate_ref_atmo:
        ref_atmo = np.zeros((num,), dtype=np.float32)
        a = np.zeros((num_time, 3), dtype=np.float32)
        a[:, 2] = 4 * np.pi / obj.wavelength  # atmospheric delay at reference acquisition
    else:
        a = np.zeros((num_time, 2))

    a[:, 1] = 4 * np.pi / obj.wavelength * tbase  # velocity

    for p in range(obj.num_points):
        coherent_time = ~np.isnan(obj.phase)[p, :]  # identify incoherent images/ifgs
        obv_vec = obj.phase[p, coherent_time]
        a[:, 0] = 4 * np.pi / obj.wavelength * pbase / (obj.slant_range[p] * np.sin(obj.loc_inc[p]))  # demerr

        x_hat = np.linalg.lstsq(a[coherent_time, :], obv_vec, rcond=None)[0]
        demerr[p] = x_hat[0]
        vel[p] = x_hat[1]
        if estimate_ref_atmo:
            ref_atmo[p] = x_hat[2]
        v_hat[p, coherent_time] = obv_vec - np.matmul(a[coherent_time, :], x_hat)
        omega[p] = np.dot(v_hat[p, coherent_time], v_hat[p, coherent_time].transpose())
        coherence[p] = np.abs(np.mean(np.exp(1j * v_hat[p, :])))

    if not estimate_ref_atmo:
        ref_atmo = None

    return vel, demerr, ref_atmo, coherence, omega, v_hat


def splitImageIntoBoxesRngAz(*, length: int, width: int, num_box_az: int, num_box_rng: int):
    """Split the image into several boxes.

    (adapted from mintpy.ifgram_inversion.split2boxes)

    Parameters
    ----------
    num_box_rng: int
        Number of boxes in range direction
    num_box_az:
        Number of boxes in azimuth direction
    length: int
        length of the image
    width: int
        width of the image

    Returns
    -------
    box_list: list
        of tuple of 4 int (xmin, ymin, xmax, ymax)
    num_box: int
        number of boxes
    """
    y_step = int(np.rint((length / num_box_rng) / 10) * 10)
    x_step = int(np.rint((width / num_box_az) / 10) * 10)

    box_list = []
    y0 = 0
    y1 = 0
    while y1 != length:
        x0 = 0
        x1 = 0
        # y1 = min([length, y0 + y_step])
        if y0 + y_step + int(np.rint(y_step / 2)) > length:
            y1 = length
        else:
            y1 = y0 + y_step
        while x1 != width:
            if x0 + x_step + int(np.rint(x_step / 2)) > width:
                x1 = width
            else:
                x1 = x0 + x_step
            # x1 = min([width, x0 + x_step])
            box = (x0, y0, x1, y1)
            box_list.append(box)
            x0 = x1
        y0 = y1

    num_box = len(box_list)
    return box_list, num_box


def preparePatches(*, num_patches: int, width: int, length: int, logger: Logger):
    """Create patches to subset the image stack for parallel processing to reduce memory usage.

    Parameters
    ----------
    num_patches: int
        number of patches to split the image into
    width: int
        width of the image
    length: int
        length of the image
    logger: Logger
        logging handler

    Returns
    -------
    box_list: list
        tuples with the radar coordinates of the boxes
    num_patches: int
        number of actual patches created by the function
    """
    patch_size_lut = {
        1: (1, 1),
        2: (1, 2),
        3: (1, 3),
        4: (2, 2),
        6: (2, 3),
        8: (2, 4),
        10: (2, 5),
        12: (3, 4),
        15: (3, 5),
        28: (4, 7),
    }
    if num_patches == 1:
        box_list = [tuple(i for i in (0, 0, width, length))]
        num_patches = 1
    else:
        num_patches = num_patches
        if num_patches > max(patch_size_lut.keys()):
            num_patches = max(patch_size_lut.keys())
            logger.info(msg=f"Number of patches is higher than expected. Reduce to {num_patches} boxes.")
        else:
            while not (num_patches in patch_size_lut.keys()):
                num_patches += 1
        box_list, num_patches = splitImageIntoBoxesRngAz(length=length,
                                                         width=width,
                                                         num_box_az=patch_size_lut[num_patches][1],
                                                         num_box_rng=patch_size_lut[num_patches][0])
        logger.info(msg=f"Process {num_patches} patches " +
                    f"({patch_size_lut[num_patches][1]} x {patch_size_lut[num_patches][0]}).")
    return box_list, num_patches


def splitDatasetForParallelProcessing(*, num_samples: int, num_cores: int):
    """Split the dataset into chunks of similar size for processing them in parallel.

    Parameters
    ----------
    num_samples: int
        number of samples to be split
    num_cores: int
        number of cores to split among

    Returns
    -------
    idx: list
        list of sample ranges for each core
    """
    rest = np.mod(num_samples, num_cores)
    avg_num_samples_per_core = int((num_samples - rest) / num_cores)
    num_samples_per_core = np.zeros((num_cores,), dtype=np.int64)
    num_samples_per_core[:] = avg_num_samples_per_core
    c = rest
    i = 0
    while c > 0:
        num_samples_per_core[i] += 1
        i += 1
        c -= 1

    idx = list()
    cur_idx = 0
    for i in range(num_cores):
        idx.append([cur_idx, cur_idx + num_samples_per_core[i]])
        cur_idx += num_samples_per_core[i]

    idx = [np.arange(s, e) for s, e in idx]
    return idx


def createSpatialGrid(*, coord_map_img: np.ndarray, length: int, width: int, grid_size: int, logger: Logger):
    """Create a spatial grid over the image.

    Parameters
    ----------
    coord_map_img: np.ndarray
        map coordinates of all image pixels
    length: int
        number of pixels in length of the image
    width: int
        number of pixels in width of the image
    grid_size: int
        size of the grid in [m]
    logger: Logger
        logging handler

    Returns
    -------
    box_list: list
        of tuples with the radar coordinates of the boxes
    num_box: int
        actual number of boxes created by the function
    """
    p0 = coord_map_img[:, 0, 0]
    p1 = coord_map_img[:, 0, -1]
    p2 = coord_map_img[:, -1, 0]
    dist_width = np.linalg.norm(p0 - p1)
    dist_length = np.linalg.norm(p0 - p2)

    if (dist_width < grid_size) or (dist_length < grid_size):
        logger.warning(f"The selected grid size ({grid_size}m) is larger than the spatial coverage of the image "
                       f"({dist_width}m x {dist_length}m).")
    num_box_az = int(np.round(dist_width / grid_size))
    num_box_rng = int(np.round(dist_length / grid_size))

    # avoid division by zero
    num_box_az = 1 if num_box_az == 0 else num_box_az
    num_box_rng = 1 if num_box_rng == 0 else num_box_rng

    # split image into different parts
    box_list, num_box = splitImageIntoBoxesRngAz(length=length, width=width,
                                                 num_box_az=num_box_az, num_box_rng=num_box_rng)
    logger.info(msg=f"{num_box} grid cells created.")
    return box_list, num_box


def selectBestPointsInGrid(*, box_list: list, quality: np.ndarray, sel_min: bool = True):
    """Select the best point inside a grid.

    If several pixel fulfill the criteria, the first one is selected.

    Parameters
    ----------
    box_list: list
        of tuples with the radar coordinates of the boxes
    quality: np.ndarray
        quality of the pixels
    sel_min: bool
        set to True if the minimum value shall be selected (default: True)

    Returns
    -------
    cand_mask_sparse: np.ndarray
        boolean mask of the selected pixels
    """
    cand_mask_sparse = np.zeros_like(quality).astype(np.bool_)

    for box in box_list:
        qual_box = quality[box[1]:box[3], box[0]:box[2]]
        if sel_min:
            idx_box = np.where(np.min(qual_box) == qual_box)
            if np.min(qual_box) == np.inf:  # no mininum value exists in this box
                continue
        else:  # max
            idx_box = np.where(np.max(qual_box) == qual_box)

        if idx_box[0].shape[0] > 1:  # more than one index might be found, due to quality(PS) = 1 in MiaplPy
            idx_box_tmp = [idx_box[0][0], idx_box[1][0]]
            idx_box = idx_box_tmp
        idx_img = (idx_box[0] + box[1], idx_box[1] + box[0])
        cand_mask_sparse[idx_img] = True
    return cand_mask_sparse


def spatiotemporalConsistency(*, coord_map: np.ndarray, phase: np.ndarray, wavelength: float, min_dist: int = 15,
                              max_dist: float = np.inf, knn: int = 50):
    """Spatiotemporal consistency proposed by Hanssen et al. (2008) and implemented in DePSI (van Leijen, 2014).

    Parameters
    ----------
    coord_map: np.ndarray
        map coordinates of the points
    phase: np.ndarray
        phase time series of the points
    wavelength: float
        radar wavelength in [m]
    min_dist: int
        minimum distance to other points in [m] (default: 15)
    max_dist: float
        maximum distance to other points in [m] (default: np.inf)
    knn: int
        number of nearest neighbors to consider (default: 50)

    Returns
    -------
    stc: np.ndarray
        spatiotemporal consistency of the points
    """
    from scipy.spatial import KDTree

    num_samples, num_time = phase.shape
    tree = KDTree(data=coord_map)

    stc = np.zeros((num_samples,), np.float64)

    for p in range(num_samples):
        dist, idx = tree.query([coord_map[p, 0], coord_map[p, 1]], k=knn)
        mask = (dist < max_dist) & (dist > min_dist) & (dist != 0)
        rho = list()
        for i in idx[mask]:
            diff = (phase[i, :-1] - phase[p, :-1]) - (phase[i, 1:] - phase[p, 1:])
            rho.append(wavelength / (4 * np.pi) * np.sqrt((1 / (num_time - 1) * np.sum(diff ** 2))))
        if not rho:
            stc[p] = np.nan
        else:
            stc[p] = np.min(rho)
    return stc


def temporalAutoCorrelation(*, residuals: np.ndarray, lag: int):
    """Compute the temporal autocorrelation for given time lag from the residuals.

    Parameters
    ----------
    residuals: np.ndarray
        residual phase time series (dim: num_points x num_time_steps)
    lag: int
        time lag used for computing the correlation

    Returns
    -------
    auto_corr: np.ndarray
        auto-correlation of each point (dim: num_points x lag)
    """
    num_points = residuals.shape[0]
    auto_corr = np.zeros((num_points, lag))
    for lag_num in range(1, lag + 1):
        for p in range(num_points):
            auto_corr[p, lag_num - 1] = abs(np.corrcoef(
                np.array([residuals[p, :-lag_num], residuals[p, lag_num:]]))[0][1])
    return auto_corr


def readPhasePatchwise(*, stack_obj: BaseStack, dataset_name: str, num_patches: int, cand_mask: np.ndarray,
                       point_id_img: np.ndarray,
                       logger: Logger):
    """Read the phase from a file in a patchwise manner to reduce memory usage.

    Parameters
    ----------
    stack_obj: BaseStack
        instance of class BaseStack
    dataset_name: str
        name of the dataset to read (e.g. 'ifgs' or 'phase')
    num_patches: int
        number of patches to split the image into
    cand_mask: np.ndarray
        boolean mask of the selected pixels
    point_id_img: np.ndarray
        image with point IDs for each pixel
    logger: Logger
        logging handler

    Returns
    -------
    phase_points: np.ndarray
        phase time series of the selected pixels
    """
    if dataset_name == "ifgs":
        length, width, num_images = stack_obj.getShape(dataset_name=dataset_name)
    elif dataset_name == "phase":  # result from miaplpy
        num_images, length, width = stack_obj.getShape(dataset_name=dataset_name)
    else:
        logger.error(f"Reading '{dataset_name}' is not supported.")
        raise NotImplementedError

    if num_patches == 1:
        phase_img = stack_obj.read(dataset_name=dataset_name)
        if dataset_name == "phase":  # result from miaplpy
            phase_img = np.moveaxis(phase_img, 0, -1)
            phase_points = phase_img[cand_mask, :]
        else:
            phase_points = np.angle(phase_img[cand_mask, :])
    else:
        box_list, num_patches = preparePatches(num_patches=num_patches,
                                               width=width,
                                               length=length,
                                               logger=logger)
        num_points = cand_mask[cand_mask].shape[0]
        phase_points = np.zeros((num_points, num_images), dtype=np.float32)
        start_idx = 0
        point_id_order = list()
        for idx in range(num_patches):
            bbox = box_list[idx]
            if dataset_name == "phase":  # result from miaplpy
                # slcStack has different order: starts with num_images. Adjust bbox (x0, y0, z0, x1, y1, z1)
                # read whole slcStack and subset to time span outside this function.
                box = (bbox[1], 0, bbox[0], bbox[3], num_images, bbox[2])
                phase_img = stack_obj.read(dataset_name=dataset_name, box=box, print_msg=False)
                phase_img = np.moveaxis(phase_img, 0, -1)
            else:
                phase_img = stack_obj.read(dataset_name=dataset_name, box=bbox, print_msg=False)
            cur_cand_mask = cand_mask[bbox[1]:bbox[3], bbox[0]:bbox[2]]

            # extract the wrapped phase for the selected pixels in the patch
            cur_num_points = cur_cand_mask[cur_cand_mask].shape[0]
            stop_idx = start_idx + cur_num_points
            if dataset_name == "phase":
                phase_points[start_idx:stop_idx, :] = phase_img[cur_cand_mask, :]  # miaplpy results are phases
            else:
                phase_points[start_idx:stop_idx, :] = np.angle(phase_img[cur_cand_mask, :])
            start_idx = stop_idx

            # store order of IDs to sort the points after loading all ifgs
            cur_point_id = point_id_img[bbox[1]:bbox[3], bbox[0]:bbox[2]]
            cur_point_id = cur_point_id[cur_cand_mask]
            point_id_order.append(cur_point_id)
            logger.info(msg="\r\033[KPatches read:\t {}/{}".format(idx + 1, num_patches))
        # reorder points to fit to the same structure for all datasets
        idx = np.argsort(np.hstack(point_id_order))
        phase_points = phase_points[idx, :]

    return phase_points


def detectValidAreas(*, bmap_obj: AmplitudeImage, logger: Logger):
    """Detect valid areas based on amplitude image.

    Parameters
    ----------
    bmap_obj: AmplitudeImage
        instance of class AmplitudeImage
    logger: Logger
        logging handler

    Returns
    -------
    mask_valid_area: np.ndarray
        boolean mask of the valid areas
    """
    bmap_obj.open()
    mask_valid_area = (10 ** (bmap_obj.background_map / 10)) > 0
    num_invalid = mask_valid_area[~mask_valid_area].shape[0]
    if num_invalid > 0:
        logger.info(msg=f"Number of invalid pixels found in image: {num_invalid}")
    return mask_valid_area


def setReferenceToPeakOfHistogram(*, phase: np.ndarray, vel: np.ndarray, num_bins: int = 100):
    """Set reference phase value to peak of the velocity histogram.

    It assumes that no velocity (i.e. stable area) is occuring most frequently.

    Parameters
    ----------
    phase: np.ndarray
        phase time series of the points
    vel: np.ndarray
        velocity of the points
    num_bins: int
        number of bins for the histogram (default: 100)

    Returns
    -------
    phase: np.ndarray
        phase time series adjusted by the new reference phase
    """
    # for TCS phase unwrapping: use only CCS points for reference phase
    mask_ccs = np.isnan(phase).sum(axis=1) == 0
    phase_ccs = phase[mask_ccs, :]
    vel = vel[mask_ccs]

    if phase_ccs.shape[0] < 40:  # the method will not give meaningfull results if too few points are available
        num_bins = 10

    # find most frequent velocity
    hist, bin_edges = np.histogram(vel, bins=num_bins, density=True)
    max_idx = np.argmax(hist)

    # find a set of points which have the most frequent velocity
    mask = (vel >= bin_edges[max_idx]) & (vel < bin_edges[max_idx + 1])

    # determine reference phase from mean of the phase time series of the selected points
    ref_phase = np.nanmean(phase_ccs[mask, :], axis=0)

    # adjust the phases by the reference sarvey
    phase -= ref_phase

    return phase


def checkIfRequiredFilesExist(*, path_to_files: str, required_files: list, logger: Logger):
    """
    Check if all required files exist from previous processing steps.

    Parameters
    ----------
    path_to_files: str
        path to the files
    required_files: list
        list of required files which are all checked
    logger: Logger
        logging handler

    Raises
    ------
    FileNotFoundError
        if a required file is missing
    """
    # loop over all required files and check if they exist, if not: raise error
    for file in required_files:
        if not exists(join(path_to_files, file)):
            logger.error(f"File from previous step(s) is missing: {file}.")
            raise FileNotFoundError
