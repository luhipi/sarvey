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

"""MTI module for SARvey."""

import argparse
import json
import os
import shutil
from json import JSONDecodeError
from os.path import join
import matplotlib
import sys
import logging
import time
from logging import Logger
from pydantic.schema import schema

from sarvey.console import printStep, printCurrentConfig, showLogoSARvey
from sarvey.processing import Processing
from sarvey.config import Config
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

    steps = range(args.start, args.stop + 1)

    config_default_dict = generateTemplateFromConfigModel()

    proc_obj = Processing(path=config.data_directories.path_outputs, config=config, logger=logger)

    printCurrentConfig(config_section=config.processing.dict(),
                       config_section_default=config_default_dict["processing"],
                       logger=logger)

    if config.phase_linking.phase_linking:
        printCurrentConfig(config_section=config.phase_linking.dict(),
                           config_section_default=config_default_dict["phase_linking"],
                           logger=logger)

    if 0 in steps:
        printStep(step=0, step_dict=STEP_DICT, logger=logger)
        printCurrentConfig(config_section=config.preparation.dict(),
                           config_section_default=config_default_dict["preparation"],
                           logger=logger)
        proc_obj.runPreparation()
    required_files = ["background_map.h5", "coordinates_utm.h5", "ifg_network.h5", "ifg_stack.h5",
                      "temporal_coherence.h5"]

    if 1 in steps:
        checkIfRequiredFilesExist(
            path_to_files=config.data_directories.path_outputs,
            required_files=required_files,
            logger=logger
        )
        printStep(step=1, step_dict=STEP_DICT, logger=logger)
        printCurrentConfig(config_section=config.consistency_check.dict(),
                           config_section_default=config_default_dict["consistency_check"],
                           logger=logger)
        proc_obj.runConsistencyCheck()
    required_files.append(["point_network.h5", "point_network_parameter.h5", "p1_ifg_wr.h5"])

    if 2 in steps:
        checkIfRequiredFilesExist(
            path_to_files=config.data_directories.path_outputs,
            required_files=required_files,
            logger=logger
        )
        printStep(step=2, step_dict=STEP_DICT, logger=logger)
        printCurrentConfig(config_section=config.unwrapping.dict(),
                           config_section_default=config_default_dict["unwrapping"],
                           logger=logger)
        if proc_obj.config.processing.temporal_unwrapping:
            proc_obj.runUnwrappingTimeAndSpace()
        else:
            proc_obj.runUnwrappingSpace()
    required_files.append(["p1_ifg_unw.h5", "p1_ts.h5"])

    if 3 in steps:
        checkIfRequiredFilesExist(
            path_to_files=config.data_directories.path_outputs,
            required_files=required_files,
            logger=logger
        )
        printStep(step=3, step_dict=STEP_DICT, logger=logger)
        printCurrentConfig(config_section=config.filtering.dict(),
                           config_section_default=config_default_dict["filtering"],
                           logger=logger)
        proc_obj.runFiltering()
    coh_value = int(config.filtering.coherence_p2 * 100)
    required_files.append(["p1_aps.h5", f"coh{coh_value}_ifg_wr.h5", f"coh{coh_value}_aps.h5"])

    if 4 in steps:
        checkIfRequiredFilesExist(
            path_to_files=config.data_directories.path_outputs,
            required_files=required_files,
            logger=logger
        )
        printStep(step=4, step_dict=STEP_DICT, logger=logger)
        printCurrentConfig(config_section=config.densification.dict(),
                           config_section_default=config_default_dict["densification"],
                           logger=logger)
        if proc_obj.config.processing.temporal_unwrapping:
            proc_obj.runDensificationTimeAndSpace()
        else:
            proc_obj.runDensificationSpace()

    logger.info(msg="SARvey MTI finished normally.")
    # close log-file to avoid problems with deleting the files
    if logger.hasHandlers():
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)
            handler.flush()
            handler.close()


def generateTemplateFromConfigModel():
    """GenerateTemplateFromConfigModel."""
    top_level_schema = schema([Config])
    top_level_dict = dict()
    for sec_name, sec_def in top_level_schema['definitions'].items():
        if sec_name == "Config":
            # substitute the class names of subsections in top_level_dict by the name of the sections in class Config
            for subsec_name, subsec_def in sec_def["properties"].items():
                top_level_dict[subsec_name] = top_level_dict.pop(subsec_def["title"])
            continue  # don't add "Config" to top_level_dict
        sec_dict = dict()
        for subsec_name, subsec_def in sec_def["properties"].items():
            if "default" not in subsec_def:
                sec_dict.update({subsec_name: None})
            else:
                sec_dict.update({subsec_name: subsec_def["default"]})
        top_level_dict.update({sec_name: sec_dict})

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
            f.write(json.dumps(default_config_dict, indent=4))
        return 0

    if args.print_config_explanation:
        top_level_schema = schema([Config])
        print(json.dumps(top_level_schema, indent=2))
        return 0

    if args.stop < args.start:
        logger.error(msg="Choose Start <= Stop!")
        raise ValueError

    if args.workdir is None:
        args.workdir = os.path.abspath(os.path.curdir)
    else:
        logger.info(msg="Working directory: {}".format(args.workdir))

    config_file_path = os.path.abspath(join(args.workdir, args.filepath))

    try:
        with open(config_file_path) as config_fp:
            config_dict = json.load(config_fp)
            config = Config(**config_dict)
    except JSONDecodeError as err:
        raise IOError(f'Failed to load the configuration json file => {err}')

    current_datetime = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
    log_filename = f"sarvey_log_{current_datetime}.log"
    logpath = config.logging.logfile_path
    if not os.path.exists(logpath):
        os.mkdir(logpath)
    file_handler = logging.FileHandler(filename=join(logpath, log_filename))
    file_handler.setFormatter(log_format)
    logger.addHandler(file_handler)

    logging_level = logging.getLevelName(config.logging.logging_level)
    logger.setLevel(logging_level)

    config.data_directories.path_outputs = os.path.abspath(join(args.workdir, config.data_directories.path_outputs))
    if config.consistency_check.spatial_mask_file_p1 is not None:
        config.consistency_check.spatial_mask_file_p1 = os.path.abspath(
            join(args.workdir, config.consistency_check.spatial_mask_file_p1))
    if config.filtering.spatial_mask_file_p2 is not None:
        config.filtering.spatial_mask_file_p2 = os.path.abspath(
            join(args.workdir, config.filtering.spatial_mask_file_p2))

    # create all necessary directories
    if not os.path.exists(config.data_directories.path_outputs):
        os.mkdir(config.data_directories.path_outputs)
    if not os.path.exists(join(config.data_directories.path_outputs, "pic")):
        os.mkdir(join(config.data_directories.path_outputs, "pic"))

    # copy config file to output directory to ensure that there is always a backup config file with latest parameters
    shutil.copy2(src=config_file_path, dst=join(config.data_directories.path_outputs, "config.json"))

    run(config=config, args=args, logger=logger)


if __name__ == '__main__':
    main()
