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
    (c) 2021 Engys Ltd

Description
    Custom find module for Ninja version. It runs a simple ninja script
    that return the ninja version. At the end it will define:
    NINJA_FOUND      - if system has Ninja
    NINJA_EXECUTABLE - path to ninja
    NINJA_VERSION    - Version string (MAJOR.MINOR.SUBMINOR)

[----------------------------------------------------------------------------]]


# since ninja if already found by emake, we kept NINJA_EXECUTABLE
# only to define NINJA_FOUND in a similar way as done for make
find_program(NINJA_EXECUTABLE ninja
  DOC "ninja command line")
mark_as_advanced(NINJA_EXECUTABLE)

set(ninjafile_content
"
rule checkVersionCommand
    command = ninja --version
build ninjaVersion.dummy: checkVersionCommand
build dummy_ninjaVersion_target: phony ninjaVersion.dummy
"
)

if(NINJA_EXECUTABLE)
    # Test Ninja with a simple build
    file(WRITE
        ${CMAKE_CURRENT_BINARY_DIR}/checkNinjaVersion/build.ninja
        ${ninjafile_content}
        )
    execute_process(
        COMMAND ninja
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/checkNinjaVersion/
        OUTPUT_VARIABLE NINJA_OUTPUT
        )
    # Get Ninja version
    if(NINJA_OUTPUT)
        string(REPLACE "\n" ";" NINJA_OUTPUT "${NINJA_OUTPUT}" )
        list(GET NINJA_OUTPUT -2 NINJA_VERSION) # -2 because the last index is
        # empty: [dummy text;version;empty]
        # check if NINJA_VERSION matches the X.X.X pattern
        if(NOT NINJA_VERSION MATCHES "([0-9]+)\\.([0-9]+)\\.([0-9]+)")
            unset(NINJA_VERSION)
        endif()
    endif()
    # Remove checkNinjaVersion after test ninja build
    file(REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/checkNinjaVersion/)
endif()

# Handle the QUIETLY and REQUIRED arguments and set NINJA_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Ninja DEFAULT_MSG NINJA_EXECUTABLE NINJA_VERSION)

# Rise all warnings here
if(NINJA_FOUND)
else()
    # Ninja did not ran properly
    if(NOT NINJA_VERSION)
        message(WARNING "Ninja simple build failed.\n")
    else()
        # we should never get here since ninja was already found by emake
        message(WARNING "Ninja not found.")
    endif()
endif()
