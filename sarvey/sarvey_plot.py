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

"""Plot module for SARvey."""
import argparse
import time
import os
from os.path import join, basename, dirname
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import logging
from logging import Logger
import sys
import cmcrameri as cmc

from mintpy.utils import ptime
from mintpy.utils.plot import auto_flip_direction

from sarvey.ifg_network import IfgNetwork
from sarvey.objects import Points, AmplitudeImage, BaseStack
from sarvey import console
from sarvey import viewer
from sarvey.config import loadConfiguration
import sarvey.utils as ut

try:
    matplotlib.use('QtAgg')
except ImportError as e:
    print(e)

EXAMPLE = """Example:
  sarvey_plot outputs/p2_coh60_ts.h5 -t                # plot average velocity and time series
  sarvey_plot outputs/p2_coh80_ts.h5 -m -a             # plot velocity map and DEM correction interactively
  sarvey_plot outputs/p2_coh80_ts.h5 -r -n 0 5         # plot residuals for image 0 to 5
  sarvey_plot outputs/p2_coh80_ifg_wr.h5 -p -n 0 1 -a  # plot wrapped phase of final point selection for interferogram 0
  sarvey_plot outputs/p1_ifg_wr.h5 -p -n 0 1 -a     # plot wrapped phase of the first order network
  sarvey_plot -i -a outputs/ifg_stack.h5            # interactively plot interferograms
  sarvey_plot -i outputs/ifg_stack.h5               # store interferograms as png files
  [...]
"""


def plotMap(*, obj_name: str, save_path: str, interactive: bool = False, input_path: str, logger: Logger):
    """Plot results from sarvey as map in radar coordinates.

    Plot the velocity map, DEM error, squared sum of residuals, temporal coherence and spatiotemporal consistency.

    Parameters
    ----------
    obj_name : str
        Path to the Points object file.
    save_path : str
        Path to the directory where the figures are saved.
    interactive : bool
        If True, the plots will be shown interactively.
    input_path : str
        Path to the inputs directory containing slcStack.h5 and geometryRadar.h5.
    logger : Logger
        Logger object.
    """
    console.showLogoSARvey(logger=logger, step="Plot map")

    scatter_size = 1

    point_obj = Points(file_path=obj_name, logger=logger)
    point_obj.open(input_path=input_path)

    bmap_obj = AmplitudeImage(file_path=join(dirname(obj_name), "background_map.h5"))
    vel, demerr, _, coherence, omega, v_hat = ut.estimateParameters(obj=point_obj, ifg_space=False)

    ax = bmap_obj.plot(logger=logger)
    sc = ax.scatter(point_obj.coord_xy[:, 1], point_obj.coord_xy[:, 0], c=demerr, s=scatter_size,
                    cmap=cmc.cm.cmaps["vanimo"])
    plt.colorbar(sc, label="[m]", pad=0.03, shrink=0.5)
    plt.title("DEM correction")
    plt.ylabel('Azimuth')
    plt.xlabel('Range')
    plt.tight_layout()
    plt.gcf().savefig(join(save_path, "map_dem_correction.png"), dpi=300)
    if interactive:
        plt.show()
    else:
        plt.close(plt.gcf())

    v_range = np.max(np.abs(vel * 100))

    ax = bmap_obj.plot(logger=logger)
    sc = ax.scatter(point_obj.coord_xy[:, 1], point_obj.coord_xy[:, 0], c=vel * 100, s=scatter_size,
                    cmap=cmc.cm.cmaps["roma"],
                    vmin=-v_range, vmax=v_range)
    plt.colorbar(sc, label="[cm / year]", pad=0.03, shrink=0.5)
    plt.title("Mean Velocity")
    plt.ylabel('Azimuth')
    plt.xlabel('Range')
    plt.tight_layout()
    plt.gcf().savefig(join(save_path, "map_velocity.png"), dpi=300)
    if interactive:
        plt.show()
    else:
        plt.close(plt.gcf())

    ax = bmap_obj.plot(logger=logger)
    sc = ax.scatter(point_obj.coord_xy[:, 1], point_obj.coord_xy[:, 0], c=coherence, vmin=0, vmax=1, s=scatter_size,
                    cmap=cmc.cm.cmaps["lajolla"])
    plt.colorbar(sc, label="[-]", pad=0.03, shrink=0.5)
    plt.title("Temporal coherence")
    plt.ylabel('Azimuth')
    plt.xlabel('Range')
    plt.tight_layout()
    plt.gcf().savefig(join(save_path, "map_coherence.png"), dpi=300)
    if interactive:
        plt.show()
    else:
        plt.close(plt.gcf())

    stc = ut.spatiotemporalConsistency(coord_utm=point_obj.coord_utm, phase=point_obj.phase,
                                       wavelength=point_obj.wavelength,
                                       min_dist=50, max_dist=np.inf, knn=40)

    ax = bmap_obj.plot(logger=logger)
    sc = ax.scatter(point_obj.coord_xy[:, 1], point_obj.coord_xy[:, 0], c=stc * 100, s=scatter_size,
                    cmap=cmc.cm.cmaps["lajolla"])
    plt.colorbar(sc, label="[cm]", pad=0.03, shrink=0.5)
    plt.title("Spatiotemporal consistency")
    plt.ylabel('Azimuth')
    plt.xlabel('Range')
    plt.tight_layout()
    plt.gcf().savefig(join(save_path, "map_spatiotemporal_consistency.png"), dpi=300)
    if interactive:
        plt.show()
    else:
        plt.close(plt.gcf())


def plotTS(*, obj_name: str, input_path: str, logger: Logger):
    """Plot the derived displacement time series.

    Parameters
    ----------
    obj_name : str
        Path to the Points object file.
    input_path : str
        Path to the inputs directory containing slcStack.h5 and geometryRadar.h5.
    logger : Logger
        Logger object.
    """
    console.showLogoSARvey(logger=logger, step="Plot time series")

    point_obj = Points(file_path=obj_name, logger=logger)
    point_obj.open(input_path=input_path)
    if point_obj.phase.shape[1] == point_obj.ifg_net_obj.num_ifgs:
        logger.warning(msg="File contains ifg phase and not phase time series. Cannot display.")
    else:
        viewer.TimeSeriesViewer(point_obj=point_obj, logger=logger, input_path=input_path)
        plt.show()


def plotPhase(*, obj_name: str, save_path: str, image_range: tuple, interactive: bool = False, input_path: str,
              logger: Logger):
    """Plot the phase of a Points object file in geographic coordinates (WGS84).

    Plot the phase of each interferogram or each image depending on the domain of the input file.

    Parameters
    ----------
    obj_name : str
        Path to the Points object file.
    save_path : str
        Path to the directory where the figures are saved.
    image_range : tuple
        Range of images to be plotted.
    interactive : bool
        If True, the plots will be shown interactively.
    input_path : str
        Path to the inputs directory containing slcStack.h5 and geometryRadar.h5.
    logger : Logger
        Logger object.
    """
    console.showLogoSARvey(logger=logger, step="Plot phase")

    point_obj = Points(file_path=obj_name, logger=logger)
    point_obj.open(input_path=input_path)

    if image_range is None:
        viewer.plotIfgs(phase=point_obj.phase, coord=point_obj.coord_lalo, ttl="Phase")
    else:
        viewer.plotIfgs(phase=point_obj.phase[:, image_range[0]:image_range[1]], coord=point_obj.coord_lalo,
                        ttl="Phase")

    plt.gcf().savefig(join(save_path, "{}_phase.png".format(basename(obj_name)[:-3])), dpi=300)

    if interactive:
        plt.show()


def plotResidualPhase(*, obj_name: str, save_path: str, image_range: tuple, interactive: bool = False,
                      input_path: str, logger: Logger):
    """Plot the residual phase of a Points object file in geographic coordinates (WGS84).

    The residuals are derived by substracting the phase contributions based on the estimated parameters.

    Parameters
    ----------
    obj_name : str
        Path to the Points object file.
    save_path : str
        Path to the directory where the figures are saved.
    image_range : tuple
        Range of images to be plotted.
    interactive : bool
        If True, the plots will be shown interactively.
    input_path : str
        Path to the inputs directory containing slcStack.h5 and geometryRadar.h5.
    logger : Logger
        Logger object.
    """
    console.showLogoSARvey(logger=logger, step="Plot residual phase")

    point_obj = Points(file_path=obj_name, logger=logger)
    point_obj.open(input_path=input_path)

    if point_obj.phase.shape[1] == point_obj.ifg_net_obj.num_ifgs:
        v_hat = ut.estimateParameters(obj=point_obj, ifg_space=True)[-1]
    else:
        v_hat = ut.estimateParameters(obj=point_obj, ifg_space=False)[-1]

    if image_range is None:
        viewer.plotIfgs(phase=v_hat, coord=point_obj.coord_lalo, ttl="Residual phase")
    else:
        viewer.plotIfgs(phase=v_hat[:, image_range[0]:image_range[1]], coord=point_obj.coord_lalo, ttl="Residual phase")

    plt.gcf().savefig(join(save_path, "{}_residual_phase.png".format(basename(obj_name)[:-3])), dpi=300)

    if interactive:
        plt.show()


def plotAllIfgs(*, obj_name: str, save_path: str, interactive: bool = False, logger: Logger):
    """Plot all interferograms inside the ifg_stack.h5 file.

    If interactivity is enabled, the plots are shown and figures are not saved. Otherwise, the figures are
    not shown but saved as png files.
    If the ifg_network.h5 file is available, the baselines are displayed in the title of each interferogram.

    Parameters
    ----------
    obj_name : str
        Path to the ifg_stack.h5 file.
    save_path : str
        Path to the directory where the figures are saved.
    interactive : bool
        If True, the plots will be shown interactively.
    logger : Logger
        Logger object.
    """
    console.showLogoSARvey(logger=logger, step="Plot interferograms")

    if obj_name.split("/")[-1] != "ifg_stack.h5":
        logger.warning(msg="Cannot plot ifgs from {}".format(obj_name))
        return

    ifg_stack_obj = BaseStack(file=obj_name, logger=logger)
    ifgs = ifg_stack_obj.read(dataset_name="ifgs")

    path_ifg_net = join(dirname(obj_name), "ifg_network.h5")
    if os.path.exists(path_ifg_net):
        ifg_net_obj = IfgNetwork()
        ifg_net_obj.open(path=path_ifg_net)
    else:
        logger.warning(msg="'ifg_network.h5' is not available in the same directory as 'ifg_stack.h5'. "
                           "No baseline information available.")
        ifg_net_obj = None

    num_ifgs = ifgs.shape[2]

    prog_bar = ptime.progressBar(maxValue=num_ifgs)
    start_time = time.time()
    logger.info(msg="plot and save figures of ifgs.")
    for i in range(num_ifgs):
        fig = plt.figure(figsize=(15, 5))
        ax = fig.add_subplot()
        ifg = np.angle(ifgs[:, :, i])
        ifg[ifg == 0] = np.nan
        im = plt.imshow(ifg, cmap=cmc.cm.cmaps["romaO"], interpolation='nearest', vmin=-np.pi, vmax=np.pi)
        auto_flip_direction(ifg_stack_obj.metadata, ax=ax, print_msg=False)
        ax.set_xlabel("Range")
        ax.set_ylabel("Azimuth")
        plt.colorbar(im, ax=ax, label="[rad]", pad=0.03, shrink=0.5)
        if ifg_net_obj is not None:
            if ifg_net_obj.dates is not None:
                date1 = ifg_net_obj.dates[ifg_net_obj.ifg_list[i][0]]
                date2 = ifg_net_obj.dates[ifg_net_obj.ifg_list[i][1]]
                ttl = "{date1} - {date2}\nbaselines: {tbase} days, {pbase} m".format(
                    date1=date1,
                    date2=date2,
                    tbase=int(np.round(ifg_net_obj.tbase_ifg[i] * 365.25)),
                    pbase=int(np.round(ifg_net_obj.pbase_ifg[i]))
                )
            else:
                ttl = "baselines: {tbase} days, {pbase} m".format(
                    tbase=int(np.round(ifg_net_obj.tbase_ifg[i] * 365.25)),
                    pbase=int(np.round(ifg_net_obj.pbase_ifg[i]))
                )
            plt.title(ttl)
        plt.tight_layout()
        if interactive:
            plt.show()
        else:
            fig.savefig(join(save_path, "{}_ifg".format(i)), dpi=300)
            plt.close(fig)
        prog_bar.update(value=i + 1, every=1, suffix='{}/{} ifgs'.format(i + 1, num_ifgs))
    prog_bar.close()
    m, s = divmod(time.time() - start_time, 60)
    logger.debug(msg='time used: {:02.0f} mins {:02.1f} secs.'.format(m, s))


def createParser():
    """Create_parser."""
    parser = argparse.ArgumentParser(
        description='Plot results from MTI\n\n',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=EXAMPLE)

    parser.add_argument('input_file', help='Path to the input file')

    parser.add_argument('-t', '--plot-ts', default=False, dest="plotTS", action="store_true",
                        help='Creates an interactive time series viewer.')

    parser.add_argument('-p', '--plot-phase', default=False, dest="plotPhase", action="store_true",
                        help='Plots the phase.')

    parser.add_argument('-r', '--plot-residual-phase', default=False, dest="plot_res_phase", action="store_true",
                        help='Plots the residual phase after substracting known components.')

    parser.add_argument('-m', '--plot-map', default=False, dest="plotMap", action="store_true",
                        help='Plots the velocity map and DEM correction.')

    parser.add_argument('-i', '--plot-all-ifgs', default=False, dest="plotAllIfgs", action="store_true",
                        help='Plots all ifgs.')

    parser.add_argument('-n', '--image_range', default=None, dest="image_range", nargs=2, type=int,
                        help='Reduces the number of phase images to the given range. Has no effect on -m and -t. Tuple.'
                             '(default: all images).')

    parser.add_argument('-a', '--interactive', default=False, dest="interactive", action="store_true",
                        help='Enables interactive visualisation of figures. Is always ON for plot-ts.')

    parser.add_argument('-w', '--workdir', default=None, dest="workdir",
                        help='Working directory (default: current directory).')

    return parser


def main(iargs=None):
    """Run Main."""
    parser = createParser()
    args = parser.parse_args(iargs)

    if args.workdir is None:
        args.workdir = os.getcwd()

    # initiate logger
    logging_level = logging.getLevelName('DEBUG')

    log_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)

    current_datetime = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
    log_filename = f"sarvey_plot_{current_datetime}.log"
    if not os.path.exists(os.path.join(os.getcwd(), "logfiles")):
        os.mkdir(os.path.join(os.getcwd(), "logfiles"))
    file_handler = logging.FileHandler(filename=os.path.join(os.getcwd(), "logfiles", log_filename))
    file_handler.setFormatter(log_format)
    logger.addHandler(file_handler)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(log_format)
    logger.addHandler(console_handler)
    logger.setLevel(logging_level)

    logger.info("Working directory: {}".format(args.workdir))
    args.input_file = join(args.workdir, args.input_file)

    config_file_path = os.path.abspath(join(args.workdir, dirname(args.input_file), "config.json"))

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

    folder_name = "p1" if "p1" in basename(args.input_file) else basename(args.input_file)[:8]
    folder_name = "ifgs" if "ifg_stack" in basename(args.input_file) else folder_name

    save_path = join(dirname(args.input_file), "pic", folder_name)
    if not os.path.exists(save_path):
        if not args.plotTS:  # not needed for interactive time series
            os.mkdir(save_path)

    selected = False
    if args.plotTS:
        # todo: read input_path from config file in same directory as file to be able to load height from geometryRadar
        plotTS(obj_name=args.input_file, input_path=config.general.input_path, logger=logger)
        selected = True

    if args.plotPhase:
        plotPhase(obj_name=args.input_file, save_path=save_path, image_range=args.image_range,
                  interactive=args.interactive, input_path=config.general.input_path, logger=logger)
        selected = True

    if args.plot_res_phase:
        plotResidualPhase(obj_name=args.input_file, save_path=save_path, image_range=args.image_range,
                          interactive=args.interactive, input_path=config.general.input_path, logger=logger)
        selected = True

    if args.plotMap:
        plotMap(obj_name=args.input_file, save_path=save_path, interactive=args.interactive,
                input_path=config.general.input_path, logger=logger)
        selected = True

    if args.plotAllIfgs:
        plotAllIfgs(obj_name=args.input_file, save_path=save_path, interactive=args.interactive, logger=logger)
        selected = True

    if not selected:
        logger.info(msg="No action chosen.")

    # close log-file to avoid problems with deleting the files
    if logger.hasHandlers():
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)
            handler.flush()
            handler.close()


if __name__ == '__main__':
    main()
