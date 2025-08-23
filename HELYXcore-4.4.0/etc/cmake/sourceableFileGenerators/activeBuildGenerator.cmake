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


# Reset the activeBuild
set(ACTIVE_BUILD_FILE ${HELYX_PROJECT_DIR}/platforms/activeBuild.shrc)
set(ACTIVE_BUILD_PRECISION_FILE ${HELYX_PROJECT_DIR}/platforms/activeBuild-${HELYX_PRECISION_OPTION}.shrc)
set(ACTIVE_BUILD_TARGET "\${platforms_dir}/${HELYX_OPTIONS}_${HELYX_MPI_NAME}.shrc")

# Always generate shrc files, because they're needed by emake
# Emake only uses the activeBuild.cmake build not the activeActive-<DP|SP>.cmake
if(EXISTS "${ACTIVE_BUILD_FILE}")
    file(READ "${ACTIVE_BUILD_FILE}" ACTIVE_BUILD_CONTENT)
    string(REGEX MATCH "\n\\.(.*)" OLD_ACTIVE_BUILD "${ACTIVE_BUILD_CONTENT}")

    string(REGEX REPLACE
        "(\n\\. *)\"*([^\"]+)"
        "\\1;\\2;" RESULT
        ${ACTIVE_BUILD_CONTENT}
        )
    list(GET RESULT 1 OLD_ACTIVE_BUILD)
    string(STRIP ${OLD_ACTIVE_BUILD} OLD_ACTIVE_BUILD)

    if(NOT "${ACTIVE_BUILD_TARGET}" STREQUAL "${OLD_ACTIVE_BUILD}")
        message(STATUS "Note:\n"
        "    The active build has changed.  You may wish to re-source your run-time environment.\n"
        "    activeBuild.shrc used to point to the following location:\n"
        "        ${OLD_ACTIVE_BUILD}\n"
        "    It now points to this location:\n"
        "        ${ACTIVE_BUILD_TARGET}\n"
        )
    endif()
endif()


set(output_string
"${CONFIG_FILE_HEADER}
#
# Description
#     This is an automatically-generated file that points to the latest rc file.
#     It's re-generated whenever userSettingsConfig.cmake is run.
#
#     Any changes made to this file will be lost when userSettingsConfig.cmake
#     is re-run, e.g. every time the CMakeCache is refreshed.
#
# [----------------------------------------------------------------------------]

# It's impossible to find the path of a sourced script in POSIX, so we must use
# specialised methods from the derived shells.  If the provided bash and zsh
# methods fail, then fall back to the absolute location at compile-time.
script_dir=\"\${BASH_SOURCE:-\${ZSH_NAME:+$0}}\"
[ -n \"$script_dir\" ] && \\
platforms_dir=$(\\cd $(dirname $script_dir) && \\pwd -L) || \\
platforms_dir=\"${HELYX_PROJECT_DIR}/platforms\"

# These are implementation details for emake.  They may be overwritten with more
# sensible values later.
export HELYX_SETTINGS_FILE=\"${HELYX_SETTINGS_FILE}\"

# shellcheck disable=SC1091
. \"${ACTIVE_BUILD_TARGET}\"
")

# Whatever else ends up in the config file, it must end in a newline to be POSIX
# compliant
string(CONCAT output_string "${output_string}" "\n")

if(NOT emake_config_file)
    message(FATAL_ERROR "emake_config_file invalid:  \"${emake_config_file}\"")
endif()

file(WRITE ${ACTIVE_BUILD_FILE} ${output_string})
file(WRITE ${ACTIVE_BUILD_PRECISION_FILE} ${output_string})