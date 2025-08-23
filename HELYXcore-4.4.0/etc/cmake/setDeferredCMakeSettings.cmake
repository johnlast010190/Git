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

function(apply_deferred_directory_definitions directory flags_list)
    message(STATUS "Applying system specific compile definitions on directory:\n\t ${directory}")
    set_property(DIRECTORY ${directory} APPEND PROPERTY COMPILE_DEFINITIONS "${flags_list}")
endfunction()


function(apply_deferred_directory_options directory flags_list)
    message(STATUS "Applying system specific compile options on directory:\n\t ${directory}")
    set_property(DIRECTORY ${directory} APPEND PROPERTY COMPILE_OPTIONS "${flags_list}")
endfunction()

function(check_of_deferred_directory_properties property_variable)

    # Set variable name for proper messages
    if("${DEFERRED_DIRECTORY_DEFINITIONS}" STREQUAL "${property_variable}")
        set(variable_name "DEFERRED_DIRECTORY_DEFINITIONS")
    elseif("${DEFERRED_DIRECTORY_OPTIONS}" STREQUAL "${property_variable}")
        set(variable_name "DEFERRED_DIRECTORY_OPTIONS")
    endif()

    # check if variable is formed by paired entries
    list(LENGTH property_variable size_list)

    math(EXPR even_odd "${size_list} % 2" OUTPUT_FORMAT DECIMAL)
    if(NOT ${even_odd} EQUAL 0)
        string(CONCAT s
            "${variable_name} variable (at userSettings) is not formed by paired sublists of diretory-flags."
            "The results may not be what expected. Review your inputs and reconfigure HELYX."
        )
        message(FATAL_ERROR "${s}")
    endif()

    # iterate over the list to get the sublist pairs of directory-flags:
    math(EXPR sublist_number "${size_list} / 2")
    foreach(index RANGE 0 ${sublist_number} 2)

        # create sublists
        list(SUBLIST property_variable ${index} 2 sublist)

        # then check if index 0 is a valid directory
        list(GET sublist 0 directory)
        if(NOT IS_DIRECTORY "${directory}")
            string(CONCAT s
                "${variable_name} variable (at userSettings): the listed directory path is not valid. "
                "The results may not be what expected. Review the path and reconfigure HELYX."
            )
            message(FATAL_ERROR "${s}")
        endif()

        # Get list of flags
        list(GET sublist 1 flags)

        # Check which SYSTEM_COMPILE_* variable is defined on CACHE,
        # then call the proper function to apply the flags
        if("${DEFERRED_DIRECTORY_DEFINITIONS}" STREQUAL "${property_variable}")
            apply_deferred_directory_definitions("${directory}" "${flags}")
        elseif("${DEFERRED_DIRECTORY_OPTIONS}" STREQUAL "${property_variable}")
            apply_deferred_directory_options("${directory}" "${flags}")
        endif()

    endforeach()

endfunction()

# check if SYSTEM_COMPILE_* is defined on CACHE.
# If defined, then call the function to apply the flags
if(NOT DEFERRED_DIRECTORY_DEFINITIONS STREQUAL "")
    check_of_deferred_directory_properties("${DEFERRED_DIRECTORY_DEFINITIONS}")
endif()
if(NOT DEFERRED_DIRECTORY_OPTIONS STREQUAL "")
    check_of_deferred_directory_properties("${DEFERRED_DIRECTORY_OPTIONS}")
endif()
