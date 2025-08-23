#[[---------------------------------------------------------------------------]
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : Dev
|    o     o     |  ENGYS Ltd. <http://engys.com/>
|       o        |
[-----------------------------------------------------------------------------]
License
    This file is part of HELYXcore.
    HELYXcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    HELYXcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HELYXcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with HELYXcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2019-2021 Engys Ltd.

Description
    Custom findModule for BISON.

[----------------------------------------------------------------------------]]


# Push CMAKE_MODULE_PATH to enable calling CMake's find_package() from here
set(STORED_CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}")
set(CMAKE_MODULE_PATH "")

# ps.: We need to force CMAKE_PREFIX_PATH to take only the BISON_ARCH_PATH for find_package,
# after that we can append everything to CMAKE_PREFIX_PATH.
# Need to do this because if system paths are defined for the ARCH_PATH of any other library
# (like system MPI in /usr/), that can take precedence and the BISON be found instead the one in BISON_ARCH_PATH
# If BISON_ARCH_PATH is empty, then system path will be search anyway

# Nothing special here - just respect BISON_ARCH_PATH
set(STORED_CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}")
set(CMAKE_PREFIX_PATH "${BISON_ARCH_PATH}")

set(min_bison_version 3.0.0)
find_package(BISON REQUIRED ${min_bison_version})

set(CMAKE_PREFIX_PATH "${STORED_CMAKE_PREFIX_PATH}")
list(APPEND CMAKE_PREFIX_PATH "${BISON_ARCH_PATH}")

# The bison find module version check doesn't work, so we have to implement our
# own.
if("${BISON_VERSION}" VERSION_LESS ${min_bison_version})
    message(WARNING
    "The Bison library (on which swak4Foam depends) was found, but the version was insufficient.
    Version found:  ${BISON_VERSION}
    Minimum required version:  ${min_bison_version}
    To specify an alternative path to a newer version of Bison, please set
    BISON_ARCH_PATH to point to the directory containing the bison \"lib\" (or
    \"lib64\")and \"bin\" directories.
    ")
    set(BISON_FOUND FALSE)
endif()

# Pop CMAKE_MODULE_PATH
set(CMAKE_MODULE_PATH "${STORED_CMAKE_MODULE_PATH}")
