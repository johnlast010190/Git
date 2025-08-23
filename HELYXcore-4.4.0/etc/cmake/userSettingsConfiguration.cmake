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
    Does some simple sanity checks of the user settings.

[----------------------------------------------------------------------------]]



# ============================================================================ #
# ----------------------------- Check variables ------------------------------ #
# ============================================================================ #

# Check things that are already set.  Note that there are very few hard
# requirements for a HELYX build (e.g. ThirdParty is not actually required), but
# cmake requirements (like output dirs being writable, etc...) can be checked.
# Proper HELYX_PROJECT_DIR check and warning is being done in: generateEmakePrerequisites.cmake
if(NOT IS_DIRECTORY "${HELYX_PROJECT_DIR}")
    # I'd be surprised if we get here - HELYX_PROJECT_DIR is probably used
    # somewhere before this is called...
    string(CONCAT s
        "HELYX_PROJECT_DIR evaluated as the following string, which is not a "
        "directory:\n"
        "\"${HELYX_PROJECT_DIR}\""
    )
    message(FATAL_ERROR "${s}")
elseif("" STREQUAL "${HELYX_COMPILER_NAME}")
    string(CONCAT s
        "HELYX_COMPILER_NAME is not set, and HELYX failed to set a sensible "
        "default\nPlease set HELYX_COMPILER_NAME in the settings file."
    )
    message(FATAL_ERROR "${s}")
elseif(NOT "${HELYX_PRECISION_OPTION}" MATCHES "[S|D]P")
    string(CONCAT s
        "HELYX_PRECISION_OPTION was set to \"${HELYX_PRECISION_OPTION}\"\n"
        "Allowed values are \"SP\" or \"DP\""
    )
    message(FATAL_ERROR "${s}")
elseif(NOT "${HELYX_LABEL_SIZE}" MATCHES "[32|64]")
    string(CONCAT s
        "HELYX_LABEL_SIZE was set to \"${HELYX_LABEL_SIZE}\"\n"
        "Allowed values are \"32\" or \"64\""
    )
    message(FATAL_ERROR "${s}")
elseif(NOT "${CMAKE_BUILD_TYPE}" IN_LIST CMAKE_CONFIGURATION_TYPES)
    string(CONCAT s
        "CMAKE_BUILD_TYPE was set to \"${CMAKE_BUILD_TYPE}\"\n"
        "Allowed values are \"${CMAKE_CONFIGURATION_TYPES}\" (see "
        "advancedSettings.cmake for more details)"
    )
    message(FATAL_ERROR "${s}")
elseif("" STREQUAL "${HELYX_OPTIONS}")
    message(FATAL_ERROR "HELYX_OPTIONS is empty")
endif()


# Checking if ThirdParty GCC/G++ is being used
set(third_party_gcc "FALSE")
get_filename_component(REAL_TP_DIR "${HELYX_THIRDPARTY_DIR}" REALPATH)
get_filename_component(REAL_CMAKE_CXX_COMPILER "${CMAKE_CXX_COMPILER}" REALPATH)

if(${REAL_CMAKE_CXX_COMPILER} MATCHES ${REAL_TP_DIR})
    string(CONCAT s
        "Using custom ThirdParty Gcc and G++. \n"
        "Please make sure that relevant compiler libraries and "
        "dependencies (gmp, mpc, mpfr) are listed in \"LD_LIBRARY_PATH\". \n"
        )
    message(WARNING "${s}")
    # then write a warning message in {HELYX_OPTION}.shrc
    set(third_party_gcc "TRUE")
endif()


# Checking if old environment variables (WM_* and FOAM_*) are defined
execute_process(
    COMMAND env
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/
    OUTPUT_VARIABLE ALL_ENV_VARIABLES
)

# make replacements to create a list of env variables
string(REPLACE ";" "," ALL_ENV_VARIABLES "${ALL_ENV_VARIABLES}")
string(REPLACE "\n" ";" ALL_ENV_VARIABLES "${ALL_ENV_VARIABLES}")
set(DEPRECATED_ENV_VARIABLES "")

string(REPLACE ";" "|" regex_whitelist "${ENVIRONMENT_VARIABLE_WARNING_WHITELIST}")
if("" STREQUAL "${regex_whitelist}")
    set(regex_whitelist "$^")  # Impossible regex, will never match
endif()

foreach(env_var IN LISTS ALL_ENV_VARIABLES)
    if("${env_var}" MATCHES "^WM_|^FOAM_"
        AND NOT "${env_var}" MATCHES "${regex_whitelist}"
    )
        list(APPEND DEPRECATED_ENV_VARIABLES "${env_var}")
    endif()
endforeach()

if(NOT "${DEPRECATED_ENV_VARIABLES}" STREQUAL "")
    string(REPLACE ";" "\n\t" DEPRECATED_ENV_VARIABLES "${DEPRECATED_ENV_VARIABLES}")
    string(CONCAT s
        "Environment variables prefixed with \"FOAM_\" or \"WM_\" were "
        "detected.\n"
        "These variables will have no effect on the HELYXcore compilation.\n"
        "Any changes to the HELYXcore compilation environment should typically "
        " be made in the user settings file at "
        "\"HELYX_PROJECT_DIR/etc/userSettings.cmake\".\n"
        "The variables detected are as follows:\n"
        "\t${DEPRECATED_ENV_VARIABLES}\n"
        )
    message(WARNING "${s}")
endif()


list(APPEND RUNTIME_SHELL_MINUS_AVAILABLE_SHELLS ${RUNTIME_SHELL})
list(REMOVE_ITEM RUNTIME_SHELL_MINUS_AVAILABLE_SHELLS sh csh python)

if(RUNTIME_SHELL_MINUS_AVAILABLE_SHELLS)
    message(WARNING
    " Unrecognised items in RUNTIME_SHELL:"
    " \"${RUNTIME_SHELL_MINUS_AVAILABLE_SHELLS}\" "
    " (recognised values are \"sh\", \"csh\")"
    )
endif()

if("OFF" STREQUAL "${SCOTCH_REQUIRED}" AND "${PTSCOTCH_REQUIRED}" MATCHES "ON|AUTO")
    message(WARNING
    "SCOTCH_REQUIRED set to \"${SCOTCH_REQUIRED}\", but PTSCOTCH_REQUIRED set to \"${PTSCOTCH_REQUIRED}\"
    Normally, scotch is required to use ptscotch.  Your configuration may not build.
    ")
endif()


# There's a a bug where HELYX_OPTIONS changes in userSettings.cmake, and a re-
# configuration is triggered independently of emake/IDE.  This can change build
# flags in the cbuild directory, and trigger complete re-builds of both the new
# and the old configurations.
if(STRICT_BINARY_DIR_CHECKING)
    get_filename_component(prev_build_dir_name ${CMAKE_BINARY_DIR} NAME)
    if(NOT "${prev_build_dir_name}" STREQUAL "${HELYX_OPTIONS}_${HELYX_MPI_NAME}")
        message(FATAL_ERROR "Mismatch between CMAKE_BINARY_DIR and HELYX_OPTIONS
        CMake has been configured in the following location:
            \"${CMAKE_BINARY_DIR}\"
        However, the expected CMAKE_BINARY_DIR is as follows:
            \"${CMAKE_SOURCE_DIR}/cbuild/${HELYX_OPTIONS}_${HELYX_MPI_NAME}\"
        Please manually re-run CMake from the correct location.

        (This error can be disabled by setting \"STRICT_BINARY_DIR_CHECKING\" to \"OFF\")")
    endif()
endif()
