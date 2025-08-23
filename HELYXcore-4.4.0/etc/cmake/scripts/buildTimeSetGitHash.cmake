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
    (c) 2019-2024 Engys Ltd.

Description
    A cmake script for configuring the global.C file with the current git hash
    at build time.  It must be called with two arguments, "inputFile" and
    "outputFile".  Additionally, configuration values must also be passed, as
    they tend otherwise to only be available at configuration time:
        cmake -DinputFile=global.C.template
            -DoutputFile=${CMAKE_CURRENT_BINARY_DIR}/global.C  \
            -DmyVar='Hello world!' \
            -P ${CMAKE_SOURCE_DIR}/etc/cmake/scripts/buildTimeSetGitHash.cmake

    The CMake function configure_file() is called with the "@ONLY" option.

[----------------------------------------------------------------------------]]


cmake_minimum_required(VERSION 3.13 FATAL_ERROR)

foreach(expected_variable inputFile outputFile)
    if(NOT ${expected_variable})
        message(SEND_ERROR
        "buildTimeSetGitHash.cmake called with incorrect parameters - ${expected_variable} evaluated to false")
    endif()
endforeach()


set(GIT_HASH "unknown")
find_package(Git)
if(Git_FOUND)
    # message(STATUS "Git found: ${GIT_EXECUTABLE}")  # find_package reports this
    if(EXISTS "${HELYX_PROJECT_DIR}/.git")
        execute_process(
            COMMAND ${GIT_EXECUTABLE} --git-dir=${HELYX_PROJECT_DIR}/.git show-ref --head HEAD 2>/dev/null
            OUTPUT_VARIABLE GIT_HASH
        )
        string(SUBSTRING ${GIT_HASH} 0 12 GIT_HASH)
    else()
        message(STATUS ".git directory not found in \$HELYX_PROJECT_DIR, unable to find git hash")
    endif()
endif()


configure_file(${inputFile} ${outputFile}
    @ONLY
    NEWLINE_STYLE UNIX
)
