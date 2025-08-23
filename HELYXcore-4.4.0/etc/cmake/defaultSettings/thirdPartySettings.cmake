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
    (c) 2019-2020 Engys Ltd.

Description
    This file sets defaults for advanced and non-critical core variables, and
    checks that certain critical variables have sensible values.

[----------------------------------------------------------------------------]]

# ---------------------------- Ancillary Methods ----------------------------- #

# Macro used mainnly to set build options in VScode and Jenkins pipelines
macro(set_thirdparty_mpi_path)

    string(CONCAT invalid_default_openmpi_version
        "THIRDPARTY_MPI_NAME set to '${THIRDPARTY_MPI_NAME}' which is not supported in ThirdParty compilation. "
        "OpenMPI-4.0.2 is build by default.\n"
        "To compile and use OpenMPI-4.1.7:\n"
        "\t- set THIRDPARTY_MPI_NAME to 'openmpi-4.1.7' in the command line;\n"
        "\t- or, in the userSettings.cmake, set MPI_ARCH_PATH to '\${THIRDPARTY_PLATFORMS_PATH_SHORT}/openmpi-4.1.7'"
        "and HELYX_MPI_NAME to 'openmpi-4.1.7'.\n"
        "ATTENTION: HELYXcore compilation against OpenMPI-4.1.7 is still in the testing fase.\n"
        #"\tTo build against a 'system' MPI, set MPI_ARCH_PATH and HELYX_MPI_NAME in the userSettings.cmake file.\n"
        )

    # Ancillary variable to set the default thirdParty_MPI path via command line
    if("" STREQUAL "${THIRDPARTY_MPI_NAME}" OR
       "openmpi-4.0.2" STREQUAL "${THIRDPARTY_MPI_NAME}")
        set(default_mpi_path "${THIRDPARTY_PLATFORMS_PATH_SHORT}/openmpi-4.0.2")
    elseif("openmpi-4.1.7" STREQUAL "${THIRDPARTY_MPI_NAME}")
        set(default_mpi_path "${THIRDPARTY_PLATFORMS_PATH_SHORT}/openmpi-4.1.7")
    else()
        message(FATAL_ERROR "${invalid_default_openmpi_version}")
    endif()

endmacro()


# ============================================================================ #
# ----------------------------- ThirdParty Paths ----------------------------- #
# ============================================================================ #

# IMPORTANT:
# Do not edit THIRDPARTY_PLATFORMS_PATH_SHORT and THIRDPARTY_PLATFORMS_PATH_LONG

# The variables used to define THIRDPARTY_PLATFORMS_PATH_*
# are set in the advancedSettings.cmake
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


# ============================================================================ #
# --------------------------- ThirdParty Settings ---------------------------- #
# ============================================================================ #

# ----------------------------------- MPI ------------------------------------ #
# Set default ThirdParty OpenMPI path
set_thirdparty_mpi_path()

# MPI is always compiled
set(MPI_ARCH_PATH
    "${default_mpi_path}"
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
