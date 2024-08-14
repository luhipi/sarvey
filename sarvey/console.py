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

"""Console module for SARvey."""
from sarvey import version
from logging import Logger


def printStep(*, step: int, step_dict: dict, logger: Logger):
    """Print the current step to console.

    Parameters
    ----------
    step: int
        current step number
    step_dict: dict
        dictionary with step numbers and names
    logger: Logger
        Logging handler
    """
    logger.info("    ---------------------------------------------------------------------------------")
    logger.info(f"                         STEP {step}:     {step_dict[step]}")
    logger.info("    ---------------------------------------------------------------------------------")


def printCurrentConfig(*, config_section: dict, config_section_default: dict, logger: Logger):
    """Print the current parameters and their default values from the config file to console.

    Parameters
    ----------
    config_section: dict
        Section of the configuration class which contains the selected parameters.
    config_section_default: dict
        Config section with default values.
    logger: Logger
        Logging handler.
    """
    shift = "    "
    logger.info(shift + f"{'Parameter':>35} {'value':>15}      {'default':>10}")
    logger.info(shift + f"{'_________':>35} {'_____':>15}      {'_______':>10}")

    for key in config_section.keys():
        default = config_section_default[key]
        default = "None" if default is None else default
        default = "True" if default is True else default
        default = "False" if default is False else default

        value = config_section[key]
        value = "None" if value is None else value
        value = "True" if value is True else value
        value = "False" if value is False else value
        if default == value:
            logger.info(shift + f"{key:>35} {value:>15}      {default:>10}")
        else:
            logger.info(shift + f"{key:>35} {value:>15} <--- {default:>10}")

    logger.info("")


def showLogoSARvey(*, logger: Logger, step: str):
    """ShowLogoSARvey.

    Parameters
    ----------
    logger: Logger
        logging handler
    step: str
        Name of the step or script which is shown on the logo.
    """
    # generate_from: http://patorjk.com/software/taag/  - font: Big, style: default
    # and https://textik.com/
    logger.info(f"SARvey version: {version.__version__} - {version.__versionalias__}, {version.__versiondate__}, "
                    f"Run: {step}")
    new_logo = rf"""                                .            _____         _____
                      +------  / \  ------  / ____|  /\   |  __ \
                      |       /  /         | (___   /  \  | |__) |_   _____ _   _
                      |      /  /           \___ \ / /\ \ |  _  /\ \ / / _ \ | | |
                      |   /\\  /  /         ____) / ____ \| | \ \ \ V /  __/ |_| |
                      |  /  \\/  /         |_____/_/    \_\_|  \_\ \_/ \___|\__, |
                      | /    \  /                                            __/ |
                      | \    / /               v{version.__version__:<5} - {version.__versionalias__:<18}  |___/
                        \\  / /...             {version.__versiondate__:<20}              |
                       / \\/ /    :...                                           |
                      /  /  /         :...     {step: <20}              |
                     /  /  /              :...                                   |
                    /  /       _______        :...                      _________|
                     \/               \______     :...     ____________/         |
                      +--------------------  \________:___/  --------------------+
    """
    print(new_logo)
