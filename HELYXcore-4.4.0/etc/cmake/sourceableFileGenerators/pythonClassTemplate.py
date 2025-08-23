# ---------------------------------------------------------------------------
# |       o        |
# |    o     o     |  HELYX (R) : Open-source CFD for Enterprise
# |   o   O   o    |  Version : 4.4.0
# |    o     o     |  ENGYS Ltd. <http://engys.com/>
# |       o        |
# ---------------------------------------------------------------------------
# License
#    This file is part of HELYXcore.
#    HELYXcore is based on OpenFOAM (R) <http://www.openfoam.org/>.
#
#    HELYXcore is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    HELYXcore is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#    for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with HELYXcore.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright
#    (c) 2021 Engys Ltd.
#
# Description
#     A template for a python class that includes useful information about
#     HELYX.  This information is provided in the hope that it might be useful
#     for users with Python-based veritical integrations of HELYX.
#
# ------------------------------------------------------------------------------

# Better to use pathlib, but that would break Python 2 compatability
import os


class HelyxInformation:
    def __init__(self, HELYX_PROJECT_DIR=""):

        if HELYX_PROJECT_DIR:
            self.HELYX_PROJECT_DIR = HELYX_PROJECT_DIR
        else:
            # Assume we're in the platforms subdirectory
            self.HELYX_PROJECT_DIR = os.path.abspath(
                os.path.realpath(__file__ + "/../../")
            )

        self.HELYX_SETTINGS_FILE = "@HELYX_SETTINGS_FILE@"

        self.HELYX_OPTIONS = "@HELYX_OPTIONS@"

        self.HELYX_SRC = os.path.join(self.HELYX_PROJECT_DIR, "src")
        self.HELYX_ETC = os.path.join(self.HELYX_PROJECT_DIR, "etc")
        self.HELYX_MODULES = os.path.join(self.HELYX_PROJECT_DIR, "modules")
        self.HELYX_APPLICATIONS = os.path.join(self.HELYX_PROJECT_DIR, "applications")
        self.HELYX_SOLVERS = os.path.join(self.HELYX_APPLICATIONS, "solvers")
        self.HELYX_UTILITIES = os.path.join(self.HELYX_APPLICATIONS, "utilities")
        self.HELYX_TUTORIALS = os.path.join(self.HELYX_PROJECT_DIR, "examples")

        self.HELYX_RUNTIME_OUTPUT_DIRECTORY = @HELYX_RUNTIME_OUTPUT_DIRECTORY@
        self.HELYX_LIBRARY_OUTPUT_DIRECTORY = @HELYX_LIBRARY_OUTPUT_DIRECTORY@

        self.HELYX_CONFIG = os.path.join(self.HELYX_PROJECT_DIR, "etc/dictData")
        self.FOAM_CONFIG = os.path.join(self.HELYX_PROJECT_DIR, "etc/dictData")

        self.HELYX_THIRDPARTY_VERSION = "@HELYX_THIRDPARTY_VERSION@"
        self.HELYX_THIRDPARTY_DIR = @HELYX_THIRDPARTY_DIR_STRING@

        self.MPI_ARCH_PATH = @MPI_DIR_STRING@

        self.HELYX_BUILD_PLATFORM = "@HELYX_BUILD_PLATFORM@"
        self.HELYX_COMPILER_LIB_ARCH = "@HELYX_COMPILER_LIB_ARCH@"
        self.HELYX_COMPILER_NAME = "@HELYX_COMPILER_NAME@"
        self.HELYX_PRECISION_OPTION = "@HELYX_PRECISION_OPTION@"
        self.HELYX_LABEL_SIZE = "@HELYX_LABEL_SIZE@"
        self.HELYX_BUILD_TYPE = "@CMAKE_BUILD_TYPE@"
        self.HELYX_PROJECT_VERSION = "@HELYX_PROJECT_VERSION@"
