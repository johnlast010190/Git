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

Description
    A cmake script that passes information from CMake to emake via
    activeBuild.shrc.

[----------------------------------------------------------------------------]]

# Allow use of the IN_LIST operator
# (This is set in cmakeConfiguration.cmake, but as this file is run as a script
# this policy must also be set here)
cmake_policy(SET CMP0057 NEW)



# ============================================================================ #
# -------------------- Find or generate user settings file ------------------- #
# ============================================================================ #

# Can't use HELYX_PROJECT_DIR here, it hasn't been read in yet
get_filename_component(DEFAULT_SETTINGS_FILE
    "${CMAKE_CURRENT_LIST_DIR}/../../userSettings.cmake"
    ABSOLUTE
    )

# First set some sensible value for HELYX_SETTINGS FILE
# Prefer file in default location
# If that's a broken symlink, then error
# If that's not present then use env var (undocumented feature)
# If that's not set, then try and generate a file.
set(HELYX_SETTINGS_FILE "")
set(SETTINGS_FILE_MESSAGE "")
if(EXISTS "${DEFAULT_SETTINGS_FILE}")
    set(HELYX_SETTINGS_FILE ${DEFAULT_SETTINGS_FILE})
    set(SETTINGS_FILE_MESSAGE
    "Using the default user settings file from the following location:
    ${HELYX_SETTINGS_FILE}")
elseif(IS_SYMLINK "${DEFAULT_SETTINGS_FILE}")
    message(FATAL_ERROR "
    A broken symlink exists at the following location:
        \"${DEFAULT_SETTINGS_FILE}\"
    Unable to unambiguously determine user settings file, stopping immediately.")
elseif(EXISTS "$ENV{HELYX_SETTINGS_FILE}")
    set(HELYX_SETTINGS_FILE "$ENV{HELYX_SETTINGS_FILE}")
    set(SETTINGS_FILE_MESSAGE
    "Using the helyx settings file specified in the HELYX_SETTINGS_FILE environment variable:
    ${HELYX_SETTINGS_FILE}")
else()
    message(STATUS "Settings file not found in \$HELYX_PROJECT_DIR/etc, attempting to generate file from environment variables...")
    include(${CMAKE_CURRENT_LIST_DIR}/generateUserSettingsFile.cmake)
    set(HELYX_SETTINGS_FILE ${DEFAULT_SETTINGS_FILE})
    set(SETTINGS_FILE_MESSAGE
    "Using automatically generated settings file at the following location:
    ${HELYX_SETTINGS_FILE}")
endif()

# Having set HELYX_SETTINGS_FILE to some sensible value, check it!
# Best to give very descriptive errors here
if(NOT EXISTS "${HELYX_SETTINGS_FILE}")
    message(FATAL_ERROR "
CMake expected to find a user settings file on the following path, which is invalid:
    ${HELYX_SETTINGS_FILE}
    ")
else()
    if(IS_SYMLINK "${HELYX_SETTINGS_FILE}")
        get_filename_component(temp ${HELYX_SETTINGS_FILE} REALPATH)
        set(SETTINGS_FILE_MESSAGE "${SETTINGS_FILE_MESSAGE}
    (This is a symlink pointing to \"${temp}\")\n"
            )
    endif()
    message(STATUS ${SETTINGS_FILE_MESSAGE})
    include(${HELYX_SETTINGS_FILE})
endif()



# ============================================================================ #
# ----------------------------- Export variables ----------------------------- #
# ============================================================================ #

# Remember, this is a local CMake variable.  It's passed to emake via
# activeBUild.shrc, and emake then feeds it in to the main CMake workflow as an
# internal cache variable on the command-line.
set(emake_config_file
    "${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}_${HELYX_MPI_NAME}.shrc"
    )

# Can't generate the full sourceable file yet, because we don't have enough
# information about MPI.  However, if nothing exists we can (and should)
# generate a more minimal file that doesn't include MPI information, PATH, or
# LD_LIBRARY_PATH.  This basically just caters for "emake -g", so there's no
# need to run this unless this is called as a script.
include(${CMAKE_CURRENT_LIST_DIR}/../sourceableFileGenerators/commonSourceableFileStrings.cmake)
if("${HELYX_SETTINGS_FILE}" IS_NEWER_THAN "${emake_config_file}")
    include(${CMAKE_CURRENT_LIST_DIR}/../sourceableFileGenerators/minimumSourceableFileGenerator.cmake)
endif()

# Always update activeBuild.shrc
include(${CMAKE_CURRENT_LIST_DIR}/../sourceableFileGenerators/activeBuildGenerator.cmake)

# A more complete sourceable file (and sourceable files for alternative shells)
# is generated after MPI has been configured, because information about MPI is
# very useful.
