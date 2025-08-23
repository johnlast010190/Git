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
    Custom find module for Make version. It runs a simple make script
    that return the make version. At the end it will define:
    MAKE_FOUND      - if system has Make
    MAKE_EXECUTABLE - path to make
    MAKE_VERSION    - Version string (MAJOR.MINOR.SUBMINOR)

[----------------------------------------------------------------------------]]


find_program (MAKE_EXECUTABLE make
  DOC "make command line")
mark_as_advanced(MAKE_EXECUTABLE)

set(makefile_content
".PHONY: all
all: \; @echo $(MAKE_VERSION)
#need := ${MAKE_VERSION_MIN}
#ok := $(filter $(need),$(firstword $(sort $(MAKE_VERSION) $(need))))
"
)

if(MAKE_EXECUTABLE)

    file(WRITE
        ${CMAKE_CURRENT_BINARY_DIR}/checkMakeVersion/Makefile
        ${makefile_content}
    )
    execute_process(
        COMMAND make -s
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/checkMakeVersion/
        OUTPUT_VARIABLE MAKE_VERSION
    )
endif()

# Handle the QUIETLY and REQUIRED arguments and set MAKE_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Make DEFAULT_MSG MAKE_EXECUTABLE MAKE_VERSION)

# Rise a warning if the minimum make version is not available
if(MAKE_FOUND)
    string(CONCAT old_make_explanation
        "Make is used by CMake to generate the build system.  Some very old "
        "versions of Make can fail to resolve dependencies properly, which "
        "breaks the process of re-building HELYX.  Clean builds *should* be "
        "okay.\n"
        "The recommended fix for this problem is to install the Ninja build "
        "system, which should be available in your package manager, and run "
        "\"emake -r\".  Ninja will be automatically detected by emake if the "
        "\"ninja\" executable is available on your PATH.\n"
    )
    if(NOT MAKE_VERSION)
        message(WARNING "Make version check failed.\n${old_make_explanation}")
    elseif(MAKE_VERSION VERSION_LESS MAKE_VERSION_MIN)
        string(CONCAT s
            "The Make version may be insufficient\n"
            "\tMinimum recommended version: ${MAKE_VERSION_MIN}\n"
            "\tDetected version: ${MAKE_VERSION}\n"
            "${old_make_explanation}"
        )
        message(WARNING "${s}")
    endif()
endif()
