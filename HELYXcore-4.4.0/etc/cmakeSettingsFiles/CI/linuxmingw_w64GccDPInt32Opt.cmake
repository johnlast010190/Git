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
    This file encodes settings for cross-compiling HELYXcore for Windows.

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
    FORCE
    )
set(HELYX_THIRDPARTY_DIR
    "${PARENT_DIR}/ThirdParty-${HELYX_THIRDPARTY_VERSION}"
    CACHE PATH
    "Path to the HELYX ThirdParty directory"
    FORCE
    )

set(HELYX_COMPILER_NAME
    "Gcc"
    CACHE PATH
    "User-friendly name for CMAKE_C_COMPILER (ONLY used for setting directory names)"
    FORCE
    )



# ============================================================================ #
# ----------------------------- Compiler Settings ---------------------------- #
# ============================================================================ #

# If the executables specified by CMAKE_C_COMPILER and CMAKE_CXX_COMPILER are
# not available on the path, please specify the full paths to the executables.
# Ps.: Preferably set the compiler full path to CMAKE_*_COMPILER variables
set(CMAKE_C_COMPILER
    "x86_64-w64-mingw32-gcc"
    CACHE PATH
    "The C compiler to use (can be a path to an executable)"
    FORCE
    )
set(CMAKE_CXX_COMPILER
    "x86_64-w64-mingw32-g++-posix"
    CACHE PATH
    "The C++ compiler to use (can be a path to an executable)"
    FORCE
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

set(CMAKE_UNITY_BUILD
    OFF
    CACHE BOOL
    "Enable unity builds"
    FORCE
    )


# ============================================================================ #
# ----------------------------- Platform Settings ---------------------------- #
# ============================================================================ #

set(HELYX_BUILD_PLATFORM linuxmingw_w64 CACHE STRING "" FORCE)

# For the Windows build, HELYX_BUILD_PLATFORM contains the OS address size, so
# we don't need to include HELYX_BUILD_PLATFORM and HELYX_COMPILER_LIB_ARCH
string(CONCAT TEMP
    "${HELYX_BUILD_PLATFORM}"
    # "${HELYX_COMPILER_LIB_ARCH}"
    "${HELYX_COMPILER_NAME}"
    "${HELYX_PRECISION_OPTION}"
    "Int${HELYX_LABEL_SIZE}"
    "${CMAKE_BUILD_TYPE}"
)
set(HELYX_OPTIONS ${TEMP} CACHE STRING
    "A string used to set paths in platforms directories"
    FORCE
    )
mark_as_advanced(HELYX_OPTIONS)

set(HELYX_SYSTEM_NAME
    "MSwindows"
    CACHE PATH
    "The target platform"
    FORCE
    )

set(CMAKE_TOOLCHAIN_FILE
    ${THIS_DIR}/cmake/toolchain-mingw64.cmake
    CACHE PATH
    "The Windows toolchain file"
    FORCE
    )

# Define _MODE_T_ here to prevent mingw typedeffing mode_t to unsigned short.
# We then manually set mode_t to unsigned int in the C++ code.
set(CMAKE_CXX_FLAGS
     "-D_FILE_OFFSET_BITS=64 -D_MODE_T_"
    CACHE STRING
    "Additional CXX flags for HELYX"
    FORCE
    )

# We must set big-obj if we want to compile in debug mode
set(CMAKE_CXX_FLAGS_DEBUG
    "-D_FILE_OFFSET_BITS=64 -D_MODE_T_ --ggdb3 -O0 -fdefault-inline -DFULLDEBUG -Wa -mbig-obj"
    CACHE STRING
    "Flags used by the C++ compiler during \"Debug\" builds."
    FORCE
    )


# ============================================================================ #
# --------------------------- ThirdParty Settings ---------------------------- #
# ============================================================================ #
include(${HELYX_PROJECT_DIR}/etc/cmake/defaultSettings/thirdPartySettings.cmake)

# Note: If you want to overwrite variables set in "thirdPartySettings.cmake",
# do it here.  Don't forget to use the "FORCE" keyword.


# ----------------------------------- MPI ------------------------------------ #
# MPI is always compiled
set(MPI_ARCH_PATH
    "${HELYX_THIRDPARTY_DIR}/platforms/linuxmingw_w64Gcc/MS-MPI-8.0"
    CACHE PATH
    "Default path on which to look for the MPI library"
    FORCE
    )
get_filename_component(FRIENDLY_MPI_NAME "${MPI_ARCH_PATH}" NAME)
set(HELYX_MPI_NAME
    "${FRIENDLY_MPI_NAME}"
    CACHE STRING
    "Human-readable name for MPI, used to set paths for MPI-dependent libraries like Pstream"
    FORCE
    )
# MPI_PATH_SUFFIXES not set.  Only the  the standard subdirectories (lib,
# lib64, inc, include, etc...) of MPI_ARCH_PATH will be searched


# ---------------------------------- SCOTCH ---------------------------------- #
set(SCOTCH_REQUIRED
    ON
    CACHE STRING
    "Use third-party scotch library"
    FORCE
    )
set(SCOTCH_ARCH_PATH
    "${HELYX_THIRDPARTY_DIR}/platforms/linuxmingw_w64GccDPInt32Opt"
    CACHE PATH
    "Default path on which to look for the SCOTCH library"
    FORCE
    )
# SCOTCH_PATH_SUFFIXES not set.  Only the  the standard subdirectories (lib,
# lib64, inc, include, etc...) of SCOTCH_ARCH_PATH will be searched


# --------------------------------- PTSCOTCH --------------------------------- #
set(PTSCOTCH_REQUIRED
    ON
    CACHE STRING
    "Use third-party scotch library"
    FORCE
    )
set(PTSCOTCH_ARCH_PATH
    "${HELYX_THIRDPARTY_DIR}/platforms/linuxmingw_w64GccDPInt32Opt"
    CACHE PATH
    "Default path on which to look for the PTSCOTCH library"
    FORCE
    )
# PTSCOTCH_PATH_SUFFIXES not set.  Only the  the standard subdirectories (lib,
# lib64, inc, include, etc...) of PTSCOTCH_ARCH_PATH will be searched

# ---------------------------------- KAHIP ----------------------------------- #
set(KAHIP_REQUIRED
    OFF
    CACHE STRING
    "Use third-party KaHIP library"
    FORCE
    )
set(KAHIP_ARCH_PATH
    "${THIRDPARTY_PLATFORMS_PATH_LONG}/kahip-3.18"
    CACHE PATH
    "Default path on which to look for the KaHIP library"
    FORCE
    )

# ---------------------------------- PARHIP ----------------------------------- #
set(PARHIP_REQUIRED
    OFF
    CACHE STRING
    "Use third-party ParHIP library"
    FORCE
    )
set(PARHIP_ARCH_PATH
    "${THIRDPARTY_PLATFORMS_PATH_LONG}/kahip-3.18"
    CACHE PATH
    "Default path on which to look for the ParHIP library"
    FORCE
    )

# ----------------------------------- REGEX ---------------------------------- #
set(REGEX_REQUIRED
    ON
    CACHE STRING
    "Use third-party regex library"
    FORCE
    )

# -------------------------------- STACK_TRACE ------------------------------- #
set(STACK_TRACE_REQUIRED
    ON
    CACHE STRING
    "Use third-party stack trace library"
    FORCE
    )

# ---------------------------------- FFTW ------------------------------------ #
set(FFTW_REQUIRED
    ON
    CACHE STRING
    "Use third-party FFTW library"
    FORCE
    )
set(FFTW_ARCH_PATH
    "${HELYX_THIRDPARTY_DIR}/platforms/linuxmingw_w64Gcc/fftw-3.3.5/"
    CACHE PATH
    "Default path on which to look for the FFTW library"
    FORCE
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
    "${HELYX_PROJECT_DIR}/../../EVTK/ext/"
    CACHE PATH
    "Default path on which to look for the VTK library"
    )
set(VTK_RENDERING_BACKEND
    OpenGL2
    CACHE STRING
    "The VTK rendering backend (one of OSMesa, OpenGL2)"
    FORCE
   )


# ============================================================================ #
# ------------------------------ HELYX modules ------------------------------- #
# ============================================================================ #

# ------------------------------- Unit Tests --------------------------------- #
set(unitTests_REQUIRED
    OFF
    CACHE STRING
    "To enable Helyx Unit tests compilation"
    FORCE
    )


# ============================================================================ #
# Set sensible defaults for ThirdParty compile switches, etc...
include(${HELYX_PROJECT_DIR}/etc/cmake/parseThirdPartySettings.cmake)
