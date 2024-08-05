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

"""Coherence module for SARvey."""
import multiprocessing
import time
import numpy as np
from numba import jit
from scipy.signal import convolve2d
from logging import Logger
from miaplpy.objects.slcStack import slcStack
from sarvey.objects import BaseStack
from sarvey.utils import convertBboxToBlock


def computeIfgsAndTemporalCoherence(*, path_temp_coh: str, path_ifgs: str, path_slc: str, ifg_array: np.ndarray,
                                    time_mask: np.ndarray, wdw_size: int, num_boxes: int, box_list: list,
                                    num_cores: int, logger: Logger):
    """ComputeIfgsAndTemporalCoherence.

    Compute the interferograms and temporal coherence from the SLC stack for a given set of (spatial) patches.

    Parameters
    ----------
    path_temp_coh : str
        Path to the temporary coherence stack. The data will be stored in this file during processing.
    path_ifgs : str
        Path to the interferograms stack. The data will be stored in this file during processing.
    path_slc : str
        Path to the SLC stack. The data will be read from this file.
    ifg_array : np.ndarray
        Array containing the indices of the reference and secondary images which are used to compute the interferograms.
    time_mask : np.ndarray
        Binary mask indicating the selected images from the SLC stack.
    wdw_size : int
        Size of the filter window. Has to be odd.
    num_boxes : int
        Number of patches to enable reading and processing of larger SLC stacks.
    box_list : list
        List containing the indices of each patch.
    num_cores : int
        Number of cores for parallel processing.
    logger : Logger
        Logger object.

    Returns
    -------
    mean_amp_img : np.ndarray
        Mean amplitude image.
    """
    start_time = time.time()
    filter_kernel = np.ones((wdw_size, wdw_size), dtype=np.float64)
    filter_kernel[wdw_size // 2, wdw_size // 2] = 0

    slc_stack_obj = slcStack(path_slc)
    slc_stack_obj.open()
    temp_coh_obj = BaseStack(file=path_temp_coh, logger=logger)
    ifg_stack_obj = BaseStack(file=path_ifgs, logger=logger)

    mean_amp_img = np.zeros((slc_stack_obj.length, slc_stack_obj.width), dtype=np.float32)
    num_ifgs = ifg_array.shape[0]

    for idx in range(num_boxes):
        logger.debug(f"processing patch {idx + 1}/{num_boxes}")
        bbox = box_list[idx]
        block2d = convertBboxToBlock(bbox=bbox)

        # read slc
        slc = slc_stack_obj.read(datasetName='slc', box=bbox, print_msg=False)
        slc = slc[time_mask, :, :]
        mean_amp = np.mean(np.abs(slc), axis=0)
        mask = mean_amp == 0
        num_zeros = np.sum(mask)
        if num_zeros > 0:
            logger.debug(f"Number of zeros in patch {idx+1} of mean amplitude image: {num_zeros}.")
            # logger.debug(f"Number of zeros in patch {idx} of mean amplitude image: {num_zeros}. Replace with NaN.")
        # TODO: if possible replace zeros with NaNs to avoid log10(0) = -inf
        # mean_amp[mask] = np.nan
        mean_amp_img[bbox[1]:bbox[3], bbox[0]:bbox[2]] = np.log10(mean_amp)

        # compute ifgs
        logger.debug(f"Compute interferograms and write to file {ifg_stack_obj.file}.")
        ifgs = computeIfgs(slc=slc, ifg_array=ifg_array)
        ifg_stack_obj.writeToFileBlock(data=ifgs, dataset_name="ifgs", block=block2d, print_msg=False)
        del slc

        # filter ifgs
        avg_neighbours = np.zeros_like(ifgs)
        if num_cores == 1:
            logger.debug("Filter interferograms using 1 core for temporal coherence estimation.")
            for i in range(num_ifgs):
                avg_neighbours[:, :, i] = convolve2d(in1=ifgs[:, :, i], in2=filter_kernel, mode='same', boundary="symm")
        else:
            logger.debug(f"Filter interferograms using {num_cores} parallel cores for temporal coherence estimation.")
            pool = multiprocessing.Pool(processes=num_cores)

            args = [(
                idx,
                ifgs[:, :, idx],
                filter_kernel) for idx in range(num_ifgs)]

            results = pool.map(func=launchConvolve2d, iterable=args)

            logger.debug(f"Retrieve {len(results)} parallel processed results.")
            for j, avg_neigh in results:
                avg_neighbours[:, :, j] = avg_neigh
            del results, args, avg_neigh

        logger.debug(f"Compute temporal coherence and write to file {temp_coh_obj.file}.")
        residual_phase = np.angle(ifgs * np.conjugate(avg_neighbours))
        del ifgs, avg_neighbours
        temp_coh = np.abs(np.mean(np.exp(1j * residual_phase), axis=2))
        temp_coh_obj.writeToFileBlock(data=temp_coh, dataset_name="temp_coh", block=block2d, print_msg=False)
        logger.debug(
            (f"Number of valid temporal coherence values in patch {idx+1}: "
             f"{np.sum(np.isfinite(temp_coh))}/{np.size(temp_coh)} with "
             f"Min: {np.nanmin(temp_coh):.2}, Max: {np.nanmax(temp_coh):.2}, Mean: {np.nanmean(temp_coh):.2}"))

        del residual_phase, temp_coh
        logger.info(f"Patches processed:\t {idx+1}/{num_boxes}")

    m, s = divmod(time.time() - start_time, 60)
    logger.debug(f"\ntime used: {m:02.0f} mins {s:02.1f} secs.\n")
    return mean_amp_img


@jit(nopython=True)
def computeIfgs(*, slc: np.ndarray, ifg_array: np.ndarray):
    """ComputeIfgs.

    Parameters
    ----------
    slc : np.ndarray
        SLC stack.
    ifg_array : np.ndarray
        Array containing the indices of the reference and secondary images which are used to compute the interferograms.

    Returns
    -------
    ifgs : np.ndarray
        Interferograms.
    """
    t, length, width = slc.shape
    num_ifgs = ifg_array.shape[0]
    ifgs = np.zeros((length, width, num_ifgs), dtype=np.complex64)

    c = 0
    for i, j in ifg_array:
        ifgs[:, :, c] = slc[i, :, :] * np.conjugate(slc[j, :, :])
        c += 1
    return ifgs


def launchConvolve2d(args: tuple):
    """LaunchConvolve2d.

    Parameters
    ----------
    args : tuple
        Tuple containing the arguments for the convolution.
        Tuple contains:

        idx : int
            Index of the processed interferogram.
        ifg : np.ndarray
            Interferogram.
        filter_kernel : np.ndarray
            Filter kernel.

    Returns
    -------
    idx : int
        Index of the processed interferogram.
    avg_neighbours : np.ndarray
        Low-pass filtered phase derived as average of neighbours.
    """
    (idx, ifg, filter_kernel) = args
    avg_neighbours = convolve2d(in1=ifg, in2=filter_kernel, mode='same', boundary="symm")
    return idx, avg_neighbours
