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
    (c) 2019-2021 Engys Ltd.

[----------------------------------------------------------------------------]]

cmake_minimum_required(VERSION 3.13)

# Make cmake file inclusion convenient.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/etc/cmake")

# TODO: There are rules about which settings can or cannot be changed before `project()` is called. It
#       wouldn't surprise me if much or all of the config system really wants to be run *before* `project()`.
project(HELYXcore)

include(projectHeader)

# Make sure module is unset here, such that all thirdParty and HELYXcore targets
# get treated as core targets.  This variable is re-set appropriately in the
# modules subdir.
unset(module)

# Add this preprocessor flag everywhere
if (HELYX_ALTERNATIVE_LIBM)
     add_compile_options(-DHELYX_ALTERNATIVE_LIBM="${HELYX_ALTERNATIVE_LIBM}")
endif ()

# Configure thirdParty and modules before source so that their presence/absence
# can be detected

add_subdirectory(thirdParty)

add_subdirectory(modules)

message(TITLE "Source file configuration")

# Messaging deferred to here for cleaner output
if (HELYX_ALTERNATIVE_LIBM)
    message(STATUS "Using maths functions from alternative libm: ${HELYX_ALTERNATIVE_LIBM}")
endif ()

add_subdirectory(src)

# All executables should link OpenFOAM
link_libraries(OpenFOAM)

add_subdirectory(applications)

if(EXISTS "${HELYX_PROJECT_DIR}/examples/CMakeLists.txt")
    add_subdirectory(examples)
endif()

# Directory property COMPILE_DEFINITIONS and COMPILE_OPTIONS can only
# be accessed after directories are processed:
include(setDeferredCMakeSettings)

message(CLEAN
"[==============================================================================]

")

# -------------------------------- Unit tests -------------------------------- #

# Must be done before the sourceable files are generated, but after modules are
# detected
if(EXISTS "${HELYX_PROJECT_DIR}/unitTests/CMakeLists.txt"
        AND NOT "MSwindows" STREQUAL "${HELYX_SYSTEM_NAME}"
        AND DEFINED unitTests_REQUIRED)

    # if unitTests_REQUIRED is defined, need to make sure it is ON|AUTO|OFF
    # if unitTests_REQUIRED is not defined, then it is automatically set to AUTO
    initialise_search_variables("unitTests")

    if("${unitTests_REQUIRED}" MATCHES "^ON|^AUTO")
        set(unitTests_FOUND TRUE CACHE BOOL "" FORCE)
        get_filename_component(real_project_dir "${HELYX_PROJECT_DIR}" REALPATH)
        add_subdirectory(${real_project_dir}/unitTests)
    endif()

endif()

# -------------------------- Packing configurations -------------------------- #

include(packingConfiguration)

# ------------------------ Sourceable file generation ------------------------ #

# Must wait until MPI is configured before re-generating the sourceable file
# Must also do this after we've configured unitTests, so unitTests output dir
# can be added to PATH if necessary
# Must not be in the project footer, because that's used in external builds

set(sourceable_file_dir "${HELYX_PROJECT_DIR}/etc/cmake/sourceableFileGenerators")
include(${sourceable_file_dir}/commonSourceableFileStrings.cmake)
include(${sourceable_file_dir}/helyxOptionsFileGenerator.cmake)

# HPC module file generation
# Always generate the environment module file.
include(${sourceable_file_dir}/HPCModuleConfig.cmake)

# Always generate extra sourceable files.  They're generated in a non-obvious
# location unless the user specifically requests them.
include(${sourceable_file_dir}/cshConverter.cmake)
include(${sourceable_file_dir}/pythonClassGenerator.cmake)

# -------------------------- emake auto-completion --------------------------- #

# Must not be in the project footer, because that's used in external builds
# Don't auto-complete directory targets (they're an implementation detail, and
# there are about 600 of them)
get_property(ALL_HELYX_LIB_TARGETS_STRING GLOBAL PROPERTY ALL_HELYX_LIB_TARGETS)
get_property(ALL_HELYX_EXE_TARGETS_STRING GLOBAL PROPERTY ALL_HELYX_EXE_TARGETS)
# get_property(ALL_HELYX_DIR_TARGETS_STRING GLOBAL PROPERTY ALL_HELYX_DIR_TARGETS)

string(REPLACE ";"
    "\n" TARGETS_STRING
    "${ALL_HELYX_LIB_TARGETS_STRING}"
    "${ALL_HELYX_EXE_TARGETS_STRING}"
    # "${ALL_HELYX_DIR_TARGETS_STRING}"
    )

# Configure emake autocompletion
configure_file(${HELYX_PROJECT_DIR}/etc/cmake/templates/template-emakeAutocomplete.sh
    ${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}/.emakeAutocomplete.sh
    @ONLY
    )

# ---------------------------------------------------------------------------- #

include(projectFooter)
