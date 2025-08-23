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
    (c) 2020-2022 Engys Ltd.

[----------------------------------------------------------------------------]]



# ============================================================================ #
# ---------------------------- Directory Settings ---------------------------- #
# ============================================================================ #

set(HELYX_PROJECT_VERSION
    "dev"
    CACHE STRING
    "HELYX version"
    FORCE
    )
set(HELYX_THIRDPARTY_VERSION
    "dev"
    CACHE STRING
    "HELYX ThirdParty version"
    FORCE
    )

# Get project path - this should work if this settings file is
# anywhere under etc/
get_filename_component(THIS_DIR "${CMAKE_CURRENT_LIST_DIR}" REALPATH)
string(REGEX REPLACE "^(.*)/etc.*$" "\\1" PARENT_DIR "${THIS_DIR}")
get_filename_component(PARENT_DIR "${PARENT_DIR}" REALPATH)
set(HELYX_PROJECT_DIR
    "${PARENT_DIR}/"
    CACHE PATH
    "Path to the top-level HELYX directory"
    FORCE
    )

# This path is defined in the Docker image of the ThirdParty build.
set(HELYX_THIRDPARTY_DIR
    "/opt/engys/ThirdParty-dev"
    CACHE PATH
    "Path to the HELYX ThirdParty directory"
    FORCE
    )


# ============================================================================ #
# ----------------------------- Compiler Settings ---------------------------- #
# ============================================================================ #

# If the executables specified by CMAKE_C_COMPILER and CMAKE_CXX_COMPILER are
# not available on the path, please specify the full paths to the executables.
# Ps.: Preferably set the compiler full path to CMAKE_*_COMPILER variables
set(CMAKE_C_COMPILER
    "gcc"
    CACHE PATH
    "The C compiler to use (can be a path to an executable)"
    FORCE
    )
set(CMAKE_CXX_COMPILER
    "g++"
    CACHE PATH
    "The C++ compiler to use (can be a path to an executable)"
    FORCE
    )

set(HELYX_PRECISION_OPTION
    "DP"
    CACHE STRING
    "Float precision, either double (\"DP\") or single precision (\"SP\")"
    FORCE
    )
set(HELYX_LABEL_SIZE
    "32"
    CACHE STRING
    "Number of bits in an OpenFOAM Label"
    FORCE
    )

set(CMAKE_BUILD_TYPE
    "Opt"
    CACHE STRING
    "Build type (usually Opt, Prof, or Debug)"
    FORCE
    )



# ============================================================================ #
# ---------------------------  Advanced Settings  ---------------------------- #
# ============================================================================ #
include(${THIS_DIR}/cmake/check_helyx_project_dir.cmake)
include(${THIS_DIR}/cmake/defaultSettings/advancedSettings.cmake)

# If you want to overwrite variables set in "advancedSettings.cmake", do it
# here.  Don't forget to use the "FORCE" keyword.


# ============================================================================ #
# --------------------------- ThirdParty Settings ---------------------------- #
# ============================================================================ #

# For historical reasons, the following paths are used as default values for
# ThirdParty libraries.  These variables will not appear in the cache.
string(CONCAT THIRDPARTY_PLATFORMS_PATH_SHORT
    "${HELYX_THIRDPARTY_DIR}/platforms/"
    "${HELYX_BUILD_PLATFORM}"
    "${HELYX_COMPILER_LIB_ARCH}"
    "${HELYX_COMPILER_NAME}"
)
string(CONCAT THIRDPARTY_PLATFORMS_PATH_LONG
    "${THIRDPARTY_PLATFORMS_PATH_SHORT}"
    "${HELYX_PRECISION_OPTION}"
    "Int${HELYX_LABEL_SIZE}"
)


# ----------------------------------- MPI ------------------------------------ #
# MPI is always compiled
set(MPI_ARCH_PATH
    "${THIRDPARTY_PLATFORMS_PATH_SHORT}/openmpi-4.0.2"
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


# ---------------------------------- SCOTCH ---------------------------------- #
set(SCOTCH_REQUIRED
    ON
    CACHE STRING
    "Use third-party scotch library"
    FORCE
    )
set(SCOTCH_ARCH_PATH
    "${THIRDPARTY_PLATFORMS_PATH_LONG}/scotch_7.0.4"
    CACHE PATH
    "Default path on which to look for the SCOTCH library"
    FORCE
    )


# --------------------------------- PTSCOTCH --------------------------------- #
set(PTSCOTCH_REQUIRED
    ON
    CACHE STRING
    "Use third-party ptscotch library"
    FORCE
    )
set(PTSCOTCH_ARCH_PATH
    "${THIRDPARTY_PLATFORMS_PATH_LONG}/scotch_7.0.4"
    CACHE PATH
    "Default path on which to look for the PTSCOTCH library"
    FORCE
    )


# ---------------------------------- KAHIP ----------------------------------- #
set(KAHIP_REQUIRED
    ON
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
    ON
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


# ----------------------------------- FFTW ----------------------------------- #
set(FFTW_REQUIRED
    ON
    CACHE STRING
    "Use third-party FFTW library"
    FORCE
    )
set(FFTW_ARCH_PATH
    "${THIRDPARTY_PLATFORMS_PATH_SHORT}/fftw-3.3.4"
    CACHE PATH
    "Default path on which to look for the FFTW library"
    FORCE
    )


# ------------------------------------ TBB ----------------------------------- #
set(TBB_REQUIRED
    OFF
    CACHE STRING
    "Use third-party TBB library"
    FORCE
    )
set(TBB_ARCH_PATH
    "${THIRDPARTY_PLATFORMS_PATH_SHORT}/oneTBB-2021.12.0"
    CACHE PATH
    "Default path on which to look for the TBB library"
    FORCE
    )


# ---------------------------------- BOOST ----------------------------------- #
set(BOOST_REQUIRED
    OFF
    CACHE STRING
    "Use third-party BOOST library"
    FORCE
    )
set(BOOST_ARCH_PATH
    "${THIRDPARTY_PLATFORMS_PATH_SHORT}/boost-1.85.0"
    CACHE PATH
    "Default path on which to look for the BOOST library"
    FORCE
    )


# ------------------------------ KOKKOS-KERNELS ------------------------------ #
set(KOKKOS_REQUIRED
    OFF
    CACHE STRING
    "Use third-party kokkoslibrary"
    FORCE
    )
set(KOKKOS_ARCH_PATH
    "${THIRDPARTY_PLATFORMS_PATH_SHORT}/kokkos-4.3.00"
    CACHE PATH
    "Default path on which to look for the kokkos library"
    FORCE
    )
set(KOKKOSKERNELS_REQUIRED
    OFF
    CACHE STRING
    "Use third-party kokkos-kernels library"
    FORCE
    )
set(KOKKOSKERNELS_ARCH_PATH
    "${THIRDPARTY_PLATFORMS_PATH_SHORT}/kokkos-kernels-4.3.00"
    CACHE PATH
    "Default path on which to look for the kokkos-kernels library"
    FORCE
    )


# ---------------------------------- OPENVDB --------------------------------- #
# OpenVDB depends on TBB, BOOST and CBLOSC
# C-Blosc is used only by OpenVDB, so is can stay here
set(CBLOSC_REQUIRED
    OFF
    CACHE STRING
    "Use third-party CBLOSC_ library"
    FORCE
    )
set(CBLOSC_ARCH_PATH
    "${THIRDPARTY_PLATFORMS_PATH_SHORT}/c-blosc-1.21.5"
    CACHE PATH
    "Default path on which to look for the CBLOSC library"
    FORCE
    )


# ---------------------------------- OPENVDB --------------------------------- #
# OpenVDB optionally depends on BOOST and CBLOSC
# So set <BOOST/CBLOSC>_REQUIRED to auto to allow their detection.
# If they are available on ThirdParty, its because (probably) OpenVDB
# was compiled against them; otherwise there will be no harm in letting them here
set(OPENVDB_REQUIRED
    OFF
    CACHE STRING
    "Use third-party OPENVDB library"
    FORCE
    )
set(OPENVDB_ARCH_PATH
    "${THIRDPARTY_PLATFORMS_PATH_SHORT}/openvdb-11.0.0"
    CACHE PATH
    "Default path on which to look for the OPENVDB library"
    FORCE
    )


# -------------------------------- OPENCASCADE ------------------------------- #
# OpenCASCADE optionally depends on TBB
set(OPENCASCADE_REQUIRED
    OFF
    CACHE STRING
    "Use third-party OPENCASCADE library"
    FORCE
    )
set(OPENCASCADE_ARCH_PATH
    "${THIRDPARTY_PLATFORMS_PATH_SHORT}/opencascade-7.8.1"
    CACHE PATH
    "Default path on which to look for the OPENCASCADE library"
    FORCE
    )


# ---------------------------------- ALEMBIC --------------------------------- #
# Alembic depends on IMath
# IMath is used only by Alembic, so is can stay here
set(IMATH_REQUIRED
    OFF
    CACHE STRING
    "Use third-party IMATH library"
    FORCE
    )
set(IMATH_ARCH_PATH
    "${THIRDPARTY_PLATFORMS_PATH_SHORT}/imath-3.1.1"
    CACHE PATH
    "Default path on which to look for the IMATH library"
    FORCE
    )
set(ALEMBIC_REQUIRED
    OFF
    CACHE STRING
    "Use third-party ALEMBIC library"
    FORCE
    )
set(ALEMBIC_ARCH_PATH
    "${THIRDPARTY_PLATFORMS_PATH_SHORT}/alembic-1.8.3"
    CACHE PATH
    "Default path on which to look for the ALEMBIC library"
    FORCE
    )


# -------------------------------- preCICE ----------------------------------- #
set(PRECICE_REQUIRED
    OFF
    CACHE STRING
    "Use third-party preCICE library"
    FORCE
    )
set(PRECICE_ARCH_PATH
    "${THIRDPARTY_PLATFORMS_PATH_SHORT}/precice-2.5.0"
    CACHE PATH
    "Default path on which to look for the preCICE library"
    FORCE
    )


# -------------------------------- RTPP - VTK -------------------------------- #
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

set(runTimePostProcessing_REQUIRED
    ON
    CACHE STRING
    "To compile RTPP. It requires EVTK"
    FORCE
    )
set(VTK_ARCH_PATH
    "/opt/engys/evtk/ext/"
    CACHE PATH
    "Default path on which to look for the VTK library"
    FORCE
    )
set(VTK_RENDERING_BACKEND
    OSMesa
    CACHE STRING
    "The VTK rendering backend (one of OSMesa, OpenGL2)"
    FORCE
    )


# ------------------------------- Unit Tests --------------------------------- #
set(unitTests_REQUIRED
    OFF
    CACHE STRING
    "To enable Helyx Unit tests compilation"
    FORCE
    )


# ---------------------------------------------------------------------------- #
# Set sensible defaults for ThirdParty compile switches, etc...
include(${HELYX_PROJECT_DIR}/etc/cmake/parseThirdPartySettings.cmake)
