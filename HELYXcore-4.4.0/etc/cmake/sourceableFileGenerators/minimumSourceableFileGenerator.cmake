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


set(output_string "")
string(CONCAT output_string "${output_string}"
"${CONFIG_FILE_HEADER}
#
# Description
#     This is an automatically-generated file that exports CMake variables as
#     environment variables.  This file exports ONLY those variables which are
#     useful at compilation time, and it does NOT define a viable environment
#     in which to run HELYX.
#
#     This file is a temporary stop-gap measure, and will be replaced when
#     CMake is run properly.
#
#     These variables are exported to the shell for convenience only, and have
#     no impact on the compilation environment.
#
# [----------------------------------------------------------------------------]


${PATH_VARIABLES_STRING}

# Add directories (apart from MPI) to PATH
export PATH=\"\\
\${HELYX_PROJECT_DIR}/bin:\\
\${HELYX_RUNTIME_OUTPUT_DIRECTORY}:\\
\${HELYX_PROJECT_DIR}/wmake:\\
\$PATH\"


${HELPER_VARIABLES_STRING}

${HELPER_MPI_VARIABLES_STRING}
")

# This is a deliberately undocumented feature, implemented as a work-around for
# adding extra dependencies to extremely old Linux platforms.  Don't use this
# unless you have to - if users have system-specific requirements, they're
# generally expected to look after them themselves.
if(NOT "" STREQUAL "${HELYX_ADDITIONAL_LD_LIBRARY_PATHS}")
    string(CONCAT output_string "${output_string}"
        "\n\nexport LD_LIBRARY_PATH=\"${HELYX_ADDITIONAL_LD_LIBRARY_PATHS}:$LD_LIBRARY_PATH\"\n"
    )
endif()

# Whatever else ends up in the config file, it must end in a newline to be POSIX
# compliant
string(CONCAT output_string "${output_string}" "\n")

file(WRITE "${emake_config_file}" "${output_string}")

message(STATUS "Generated minimum source-able file at the following location:
    ${emake_config_file}
A more complete file (including MPI configuration) will be generated when the cache is refreshed.")

