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
    (c) 2021 Engys Ltd.

Description
    This is the HELYX config file.  It is a minimal file that satisfies only
    the CMake requirements for using HELYX as an imported library (e.g. target
    imports, standard CMake macros, component checking, etc...).

    Users may also require the HELYX CMake configuration.  This can be accessed
    by including the projectHeader.cmake file, which should be done separately.
[----------------------------------------------------------------------------]]

# Don't use configure_package_config_file.  We're not defining any paths in
# this file, so it adds no value.

macro(message_if_not_quiet)
    if(NOT "${HELYX_FIND_QUIETLY}")
        message(${ARGN})
    # else()
    #   No-op
    endif()
endmacro()


macro(set_and_check _var _file)
    set(${_var} "${_file}")
    if(NOT EXISTS "${_file}")
        message_if_not_quiet(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
    endif()
endmacro()


# This file should always live in platforms/${HELYX_OPTIONS}
get_filename_component(PACKAGE_PREFIX_DIR  "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
set_and_check(HELYX_PROJECT_DIR "${PACKAGE_PREFIX_DIR}")

set(HELYX-core-target-file "${CMAKE_CURRENT_LIST_DIR}/cmake/HELYX-core-targets.cmake")
if (EXISTS "${HELYX-core-target-file}")
    include("${HELYX-core-target-file}")
    message_if_not_quiet(STATUS "Using HELYX-core target file from the following location:
    \"${HELYX-core-target-file}\""
    )
else()  # Not found and not required
    set(HELYX_FOUND FALSE)
    set(HELYX_NOT_FOUND_MESSAGE
        "Could not find HELYX-core target file at this location:\n\t${HELYX-core-target-file}"
        )
endif()


# Check required components
set(_available_components @available_components@)
string(REPLACE ";"  ", " _available_components_str "${_available_components}")

foreach(comp ${HELYX_FIND_COMPONENTS})
    if(NOT "${comp}" IN_LIST _available_components)
        if(HELYX_FIND_REQUIRED_${comp})
            set(HELYX_FOUND FALSE)
            set(HELYX_NOT_FOUND_MESSAGE
            "The component \"${comp}\" is marked REQUIRED and is not available
            Available components are as follows:
                \"${_available_components_str}\""
            )
            return()
        else()  # Not found and not required
            message_if_not_quiet(STATUS
            "The component \"${comp}\" is not available
            Available components are as follows:
                \"${_available_components_str}\""
            )
        endif()
    endif()

    set(comp_target_file "${CMAKE_CURRENT_LIST_DIR}/cmake/${comp}-targets.cmake")
    # TODO:  Review this.  Suspect it uses absolute paths.
    if (EXISTS "${comp_target_file}")
        include("${comp_target_file}")
        set(HELYX_${comp}_FOUND TRUE)
    elseif(HELYX_FIND_REQUIRED_${comp})
        set(HELYX_FOUND FALSE)
        set(HELYX_NOT_FOUND_MESSAGE
            "Could not find ${comp}, the following file is required and not found:
            \"${comp_target_file}\""
            )
        return()
    else()  # Not found and not required
        message_if_not_quiet(STATUS
            "Could not find ${comp}, the following file is required and not found:
            \"${comp_target_file}\""
            )
    endif()

    if(NOT HELYX_${comp}_FOUND)
        if(HELYX_FIND_REQUIRED_${comp})
            set(HELYX_FOUND FALSE)
            set(HELYX_NOT_FOUND_MESSAGE "Unknown error finding ${comp}")
            return()
        else()  # Not found and not required
            message_if_not_quiet(STATUS "Unknown error finding ${comp}")
        endif()
    endif()
endforeach()

