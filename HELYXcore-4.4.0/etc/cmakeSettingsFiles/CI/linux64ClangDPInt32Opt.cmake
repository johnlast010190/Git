#[[---------------------------------------------------------------------------]
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
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
    (c) 2020-2024 Engys Ltd.

Description
    This file encodes settings for the clang Linux build.

    Note:  Some variables in this file are not forced.  This can lead to
    unexpected behaviours.  Use of this file for anything but automated testing
    is therefore not recommended.

[----------------------------------------------------------------------------]]


# ============================================================================ #
# ---------------------------- Directory Settings ---------------------------- #
# ============================================================================ #

set(HELYX_PROJECT_VERSION
    "4.4.0"
    CACHE STRING
    "HELYX version"
    )
set(HELYX_THIRDPARTY_VERSION
    "4.4.0"
    CACHE STRING
    "HELYX ThirdParty version"
    )

# Get project path - this should work if this settings file is
# anywhere under etc/
get_filename_component(THIS_DIR "${CMAKE_CURRENT_LIST_DIR}" REALPATH)
string(REGEX REPLACE "^(.*)/etc.*$" "\\1/.." PARENT_DIR "${THIS_DIR}")
get_filename_component(PARENT_DIR "${PARENT_DIR}" REALPATH)

set(HELYX_PROJECT_DIR
    "${PARENT_DIR}/HELYXcore-${HELYX_PROJECT_VERSION}"
    CACHE PATH
    "Path to the top-level HELYX directory"
    )
set(HELYX_THIRDPARTY_DIR
    "${PARENT_DIR}/ThirdParty-${HELYX_THIRDPARTY_VERSION}"
    CACHE PATH
    "Path to the HELYX ThirdParty directory"
    )

set(HELYX_COMPILER_NAME
    "Clang"
    CACHE STRING
    "Friendly name for the compiler, only used to set paths"
    )


# ============================================================================ #
# ----------------------------- Compiler Settings ---------------------------- #
# ============================================================================ #

# If the executables specified by CMAKE_C_COMPILER and CMAKE_CXX_COMPILER are
# not available on the path, please specify the full paths to the executables.
# Ps.: Preferably set the compiler full path to CMAKE_*_COMPILER variables
set(CMAKE_C_COMPILER
    "clang"
    CACHE PATH
    "The C compiler to use (can be a path to an executable)"
    )
set(CMAKE_CXX_COMPILER
    "clang++"
    CACHE PATH
    "The C++ compiler to use (can be a path to an executable)"
    )

set(HELYX_PRECISION_OPTION
    "DP"
    CACHE STRING
    "Float precision, either double (\"DP\") or single precision (\"SP\")"
    )
set(HELYX_LABEL_SIZE
    "32"
    CACHE STRING
    "Number of bits in an OpenFOAM Label"
    )

set(CMAKE_BUILD_TYPE
    "Opt"
    CACHE STRING
    "Build type (usually Opt, Prof, or Debug)"
    )


# ============================================================================ #
# ---------------------------  Advanced Settings  ---------------------------- #
# ============================================================================ #
include(${THIS_DIR}/cmake/check_helyx_project_dir.cmake)
include(${THIS_DIR}/cmake/defaultSettings/advancedSettings.cmake)

# Note: If you want to overwrite variables set in "advancedSettings.cmake",
# do it here.  Don't forget to use the "FORCE" keyword.


# ============================================================================ #
# --------------------------- ThirdParty Settings ---------------------------- #
# ============================================================================ #
include(${HELYX_PROJECT_DIR}/etc/cmake/defaultSettings/thirdPartySettings-unforced.cmake)

# Note: If you want to overwrite variables set in "thirdPartySettings.cmake",
# do it here.


# ------------------------------------ TBB ----------------------------------- #
set(TBB_REQUIRED
    AUTO
    CACHE STRING
    "Use third-party TBB library"
    )


# ---------------------------------- BOOST ----------------------------------- #
set(BOOST_REQUIRED
    AUTO
    CACHE STRING
    "Use third-party BOOST library"
    )


# ---------------------------------- OPENVDB --------------------------------- #
# OpenVDB depends on TBB, BOOST and CBLOSC
# C-Blosc is used only by OpenVDB, so is can stay here
set(CBLOSC_REQUIRED
    AUTO
    CACHE STRING
    "Use third-party CBLOSC_ library"
    )


# ---------------------------------- ALEMBIC --------------------------------- #
# Alembic depends on IMath
# IMath is used only by Alembic, so is can stay here
set(IMATH_REQUIRED
    AUTO
    CACHE STRING
    "Use third-party IMATH library"
    )


# ============================================================================ #
# ---------------------------- EVTK/RTPP Settings ---------------------------- #
# ============================================================================ #

# Engys fork of VTK (EVTK) is a dependency of Run-Time Post-Processing (RTPP).
# Therefore, RTPP will only compile only if EVTK is found.
# By default, RTPP compilation is disabled. To enable RTPP compilation set
# runTimePostProcessing_REQUIRED to "ON" or "AUTO". When enabled, RTPP will look
# for EVTK and if the requirements are met, RTPP will compile.
#
# runTimePostProcessing_REQUIRED not set: defaults to "AUTO";
#
# VTK_REQUIRED not set: defaults to "AUTO"; it will search for EVTK using the
# VTK_ARCH_PATH, and the HELYX GUI EVTK path ( ~/Engys/HELYX/**/GUI/ext/ );
#
# VTK_ARCH_PATH not set: only HELYX GUI EVTK path ( ~/Engys/HELYX/**/GUI/ext/ )
# will be searched and used if available
#
# VTK_PATH_SUFFIXES not set:  Only the standard subdirectories
# (lib, lib64, inc, include, etc...) of VTK_ARCH_PATH will be searched.
#
# For custom EVTK compilation, uncomment the lines below, and set the
# VTK_ARCH_PATH to the compiled EVTK

set(VTK_ARCH_PATH
    "${HELYX_PROJECT_DIR}/../EVTK/ext/"
    CACHE PATH
    "Default path on which to look for the VTK library"
    )
set(VTK_RENDERING_BACKEND
    OSMesa
    CACHE STRING
    "The VTK rendering backend (one of OSMesa, OpenGL2)"
   )


# ============================================================================ #
# ------------------------------ HELYX modules ------------------------------- #
# ============================================================================ #

# ------------------------------- Unit Tests --------------------------------- #
set(unitTests_REQUIRED
    OFF
    CACHE STRING
    "To enable Helyx Unit tests compilation"
    )


# ============================================================================ #
# Set sensible defaults for ThirdParty compile switches, etc...
include(${HELYX_PROJECT_DIR}/etc/cmake/parseThirdPartySettings.cmake)
