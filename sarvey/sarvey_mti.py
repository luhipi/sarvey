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

"""MTI module for SARvey."""

import argparse
import os
import shutil
from os.path import join

import json5
import matplotlib
import sys
import logging
import time
from logging import Logger
from pydantic import TypeAdapter

from sarvey import version
from sarvey.console import printStep, printCurrentConfig, showLogoSARvey
from sarvey.processing import Processing
from sarvey.config import Config, loadConfiguration
from sarvey.utils import checkIfRequiredFilesExist

try:
    matplotlib.use('QtAgg')
except ImportError as e:
    print(e)

EXAMPLE = """Example:
  sarvey -f config.json 0 0 -g                   # create default config file with the name config.json and exit
  sarvey -f config.json 0 0                      # run only preparation step
  sarvey -f config.json 0 4                      # run all processing steps

  sarvey -f config.json 0 0 -p                   # print explanation of the configuration parameters to console
"""

STEP_DICT = {
    0: "PREPARATION",
    1: "CONSISTENCY CHECK",
    2: "UNWRAPPING",
    3: "FILTERING",
    4: "DENSIFICATION",
}


def run(*, config: Config, args: argparse.Namespace, logger: Logger):
    """Run the specified processing steps.

    Parameters
    ----------
    config: Config
        object of configuration class.
    args: argparse.Namespace
        command line input arguments
    logger: Logger
        Logging handler.
    """
    showLogoSARvey(logger=logger, step="MTInSAR")

    start_time = time.time()

    steps = range(args.start, args.stop + 1)

    config_default_dict = generateTemplateFromConfigModel()

    proc_obj = Processing(path=config.general.output_path, config=config, logger=logger)

    printCurrentConfig(config_section=config.general.model_dump(),
                       config_section_default=config_default_dict["general"],
                       logger=logger)

    if config.phase_linking.use_phase_linking_results:
        printCurrentConfig(config_section=config.phase_linking.model_dump(),
                           config_section_default=config_default_dict["phase_linking"],
                           logger=logger)

    if 0 in steps:
        start_time_step = time.time()
        printStep(step=0, step_dict=STEP_DICT, logger=logger)
        printCurrentConfig(config_section=config.preparation.model_dump(),
                           config_section_default=config_default_dict["preparation"],
                           logger=logger)
        proc_obj.runPreparation()
        m, s = divmod(time.time() - start_time_step, 60)
        logger.info(f"Finished step 0 PREPARATION normally in {m:02.0f} mins {s:02.1f} secs.")
    required_files = ["background_map.h5", "coordinates_utm.h5", "ifg_network.h5", "ifg_stack.h5",
                      "temporal_coherence.h5"]

    if 1 in steps:
        start_time_step = time.time()
        checkIfRequiredFilesExist(
            path_to_files=config.general.output_path,
            required_files=required_files,
            logger=logger
        )
        printStep(step=1, step_dict=STEP_DICT, logger=logger)
        printCurrentConfig(config_section=config.consistency_check.model_dump(),
                           config_section_default=config_default_dict["consistency_check"],
                           logger=logger)
        proc_obj.runConsistencyCheck()
        m, s = divmod(time.time() - start_time_step, 60)
        logger.info(f"Finished step 1 CONSISTENCY CHECK normally in {m:02.0f} mins {s:02.1f} secs.")
    required_files.extend(["point_network.h5", "point_network_parameter.h5", "p1_ifg_wr.h5"])

    if 2 in steps:
        start_time_step = time.time()
        checkIfRequiredFilesExist(
            path_to_files=config.general.output_path,
            required_files=required_files,
            logger=logger
        )
        printStep(step=2, step_dict=STEP_DICT, logger=logger)
        printCurrentConfig(config_section=config.unwrapping.model_dump(),
                           config_section_default=config_default_dict["unwrapping"],
                           logger=logger)
        if proc_obj.config.general.apply_temporal_unwrapping:
            proc_obj.runUnwrappingTimeAndSpace()
        else:
            proc_obj.runUnwrappingSpace()
        m, s = divmod(time.time() - start_time_step, 60)
        logger.info(f"Finished step 2 UNWRAPPING normally in {m:02.0f} mins {s:02.1f} secs.")
        required_files.extend(["p1_ifg_unw.h5", "p1_ts.h5"])

    if 3 in steps:
        start_time_step = time.time()
        checkIfRequiredFilesExist(
            path_to_files=config.general.output_path,
            required_files=required_files,
            logger=logger
        )
        printStep(step=3, step_dict=STEP_DICT, logger=logger)
        printCurrentConfig(config_section=config.filtering.model_dump(),
                           config_section_default=config_default_dict["filtering"],
                           logger=logger)
        proc_obj.runFiltering()
        m, s = divmod(time.time() - start_time_step, 60)
        logger.info(f"Finished step 3 FILTERING normally in {m:02.0f} mins {s:02.1f} secs.")
    coh_value = int(config.filtering.coherence_p2 * 100)
    required_files.extend(["p1_aps.h5", f"p2_coh{coh_value}_ifg_wr.h5", f"p2_coh{coh_value}_aps.h5"])

    if 4 in steps:
        start_time_step = time.time()
        checkIfRequiredFilesExist(
            path_to_files=config.general.output_path,
            required_files=required_files,
            logger=logger
        )
        printStep(step=4, step_dict=STEP_DICT, logger=logger)
        printCurrentConfig(config_section=config.densification.model_dump(),
                           config_section_default=config_default_dict["densification"],
                           logger=logger)
        if proc_obj.config.general.apply_temporal_unwrapping:
            proc_obj.runDensificationTimeAndSpace()
        else:
            proc_obj.runDensificationSpace()
        m, s = divmod(time.time() - start_time_step, 60)
        logger.info(f"Finished step 4 DENSIFICATION normally in {m:02.0f} mins {s:02.1f} secs.")

    m, s = divmod(time.time() - start_time, 60)
    logger.info(f"SARvey MTI finished normally in in {m:02.0f} mins {s:02.1f} secs.")
    # close log-file to avoid problems with deleting the files
    if logger.hasHandlers():
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)
            handler.flush()
            handler.close()


def generateTemplateFromConfigModel():
    """GenerateTemplateFromConfigModel."""
    top_level_dict = {}

    for sec_name, field in Config.model_fields.items():
        sec_cls = field.annotation
        sec_dict = {}
        for subsec_name, subsec_def in sec_cls.model_fields.items():
            if subsec_def.default is not None:
                default = subsec_def.default
            else:
                default = None

            sec_dict[subsec_name] = default
        top_level_dict[sec_name] = sec_dict

    return top_level_dict


def createParser():
    """Create_parser.

    :return:
    """
    parser = argparse.ArgumentParser(
        description='Multitemporal InSAR processing workflow\n\n' +
                    'Run the following steps:\n' +
                    '0 - preparation\n' +
                    '1 - consistency check\n' +
                    '2 - spatial unwrapping\n' +
                    '3 - filtering\n' +
                    '4 - densification',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=EXAMPLE)

    parser.add_argument('start', choices={0, 1, 2, 3, 4}, type=int,
                        help='Start of processing')

    parser.add_argument('stop', choices={0, 1, 2, 3, 4}, type=int,
                        help='Stop of processing')

    parser.add_argument("-f", "--filepath", type=str, required=True, metavar="FILE",
                        help="Path to the config.json file.")

    parser.add_argument("-g", "--generate_config", action="store_true", default=False, dest="generate_config",
                        help="Write default configuration to file specified by '-f'.")

    parser.add_argument("-p", "--print_config_explanation", action="store_true", default=False,
                        dest="print_config_explanation",
                        help="Prints exhaustive explanations about configuration to console.")

    parser.add_argument('-w', '--workdir', default=None, dest="workdir",
                        help='Working directory (default: current directory).')

    parser.add_argument('--version', action='version',
                        version=f"SARvey version {version.__version__} - {version.__versionalias__}, "
                                f"{version.__versiondate__}")

    return parser


def main(iargs=None):
    """Run Main.

    :param iargs:
    """
    parser = createParser()
    args = parser.parse_args(iargs)

    # initiate logger
    logging_level = logging.getLevelName('DEBUG')  # set a default value before until level is read from config

    log_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(log_format)
    logger.addHandler(console_handler)
    logger.setLevel(logging_level)

    if args.generate_config:
        logger.info(msg=f"Write default config to file: {args.filepath}.")
        default_config_dict = generateTemplateFromConfigModel()
        with open(args.filepath, "w") as f:
            f.write(json5.dumps(default_config_dict, indent=4))
        return 0

    if args.print_config_explanation:
        top_level_schema = TypeAdapter(Config).json_schema()
        print(json5.dumps(top_level_schema, indent=2))
        return 0

    if args.stop < args.start:
        msg = f"Selected Start step ({args.start}) must be less than or equal to Stop step ({args.stop}). Exiting!"
        logger.error(msg)
        raise ValueError(msg)

    if args.workdir is None:
        args.workdir = os.path.abspath(os.path.curdir)
    logger.info(f"Working directory: {args.workdir}")

    config_file_path = os.path.abspath(join(args.workdir, args.filepath))

    config = loadConfiguration(path=config_file_path)

    current_datetime = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
    log_filename = f"sarvey_{current_datetime}.log"
    logpath = config.general.logfile_path
    if not os.path.exists(logpath):
        os.mkdir(logpath)
    file_handler = logging.FileHandler(filename=join(logpath, log_filename))
    file_handler.setFormatter(log_format)
    file_logging_level = logging.getLevelName("DEBUG")
    file_handler.setLevel(file_logging_level)
    logger.addHandler(file_handler)

    logging_level = logging.getLevelName(config.general.logging_level)
    console_handler.setLevel(logging_level)

    config.general.output_path = os.path.abspath(join(args.workdir, config.general.output_path))
    if config.consistency_check.mask_p1_file is not None:
        config.consistency_check.mask_p1_file = os.path.abspath(
            join(args.workdir, config.consistency_check.mask_p1_file))
    if config.filtering.mask_p2_file is not None:
        config.filtering.mask_p2_file = os.path.abspath(
            join(args.workdir, config.filtering.mask_p2_file))

    # create all necessary directories
    if not os.path.exists(config.general.output_path):
        os.mkdir(config.general.output_path)
    if not os.path.exists(join(config.general.output_path, "pic")):
        os.mkdir(join(config.general.output_path, "pic"))

    # copy config file to output directory to ensure that there is always a backup config file with latest parameters
    shutil.copy2(src=config_file_path, dst=join(config.general.output_path, "config.json"))

    run(config=config, args=args, logger=logger)


if __name__ == '__main__':
    main()
