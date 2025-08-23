#[[---------------------------------------------------------------------------]
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : Dev
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
    (c) 2022-2023 Engys Ltd.

Description
    Defines custom CMake functions specific to the swak4Foam module

[----------------------------------------------------------------------------]]


# ============================================================================ #
# ---------------------------------  Bison  ---------------------------------- #
# ============================================================================ #

macro(add_custom_commands_for_bison file_name)
    if(NOT ${ARGC} EQUAL 1)
        string(CONCAT s
            "add_custom_commands_for_bison() requires exactly one argument, "
            "but ${ARGC} were provided"
        )
        message(FATAL_ERROR "${s}")
    endif()

    # Make a directory to contain the generated output.
    set(WDIR "${CMAKE_CURRENT_BINARY_DIR}/${file_name}")
    file(MAKE_DIRECTORY "${WDIR}/include")

    add_custom_command(
        OUTPUT
            ${WDIR}/${file_name}.tab.cc
            ${WDIR}/include/${file_name}.tab.hh
            ${WDIR}/include/${file_name}\_location.hh
            ${WDIR}/include/${file_name}\_stack.hh
            ${WDIR}/include/${file_name}\_position.hh
        DEPENDS
            ./${file_name}.yy

        # Note: The 'Wno-deprecated' flag suppresses warnings for a deprecated
        # directive whose replacement is not supported in Bison versions < 3.3.
        COMMAND ${BISON_EXECUTABLE} -ra -v -d -Wno-deprecated ${CMAKE_CURRENT_SOURCE_DIR}/${file_name}.yy
        COMMAND sed -i -Ee "s|position.hh|${file_name}_position.hh|g" location.hh
        COMMAND cmake -E rename location.hh ${file_name}_location.hh
        COMMAND cmake -E rename stack.hh ${file_name}_stack.hh
        COMMAND cmake -E rename position.hh ${file_name}_position.hh
        COMMAND sed -i -Ee "s|stack.hh|${file_name}_stack.hh|g;s|location.hh|${file_name}_location.hh|g" ${file_name}.tab.hh

        # The cmake copy <files> ... <directory> signature was added in CMake
        # 3.5
        COMMAND cmake -E copy ${file_name}.tab.hh ${WDIR}/include/
        COMMAND cmake -E copy ${file_name}_location.hh ${WDIR}/include/
        COMMAND cmake -E copy ${file_name}_stack.hh ${WDIR}/include/
        COMMAND cmake -E copy ${file_name}_position.hh ${WDIR}/include/
        WORKING_DIRECTORY ${WDIR}
        VERBATIM
    )
    add_custom_target(bison_${file_name} DEPENDS
        ${WDIR}/${file_name}.tab.cc
        ${WDIR}/include/${file_name}.tab.hh
        ${WDIR}/include/${file_name}\_location.hh
        ${WDIR}/include/${file_name}\_stack.hh
        ${WDIR}/include/${file_name}\_position.hh
    )

    set(BISON_${file_name}_OUTPUT_SOURCE ${WDIR}/${file_name}.tab.cc)
    set(BISON_${file_name}_OUTPUT_INCLUDE_DIR ${WDIR}/include)
    set(BISON_${file_name}_OUTPUT_HEADER ${WDIR}/include/${file_name}.tab.hh)
    set(BISON_${file_name}_OUTPUTS ${WDIR}/${file_name}.tab.cc)

endmacro()

macro(target_add_bison_sources TARGET_NAME)
    foreach(source_file ${ARGN})
        add_custom_commands_for_bison("${source_file}")
        add_custom_target(${source_file}
            DEPENDS
                ${BISON_${file_name}_OUTPUT_SOURCE}
                ${BISON_${file_name}_OUTPUT_HEADER}
        )
        add_dependencies(${TARGET_NAME}_obj "${source_file}")
        add_dependencies(${TARGET_NAME}_obj bison_${source_file})

        target_sources(${TARGET_NAME}_obj
            PRIVATE "${BISON_${source_file}_OUTPUT_SOURCE}"
        )
        target_include_directories(${TARGET_NAME}_obj
            PUBLIC $<BUILD_INTERFACE:${BISON_${source_file}_OUTPUT_INCLUDE_DIR}>
        )
        target_include_directories(${TARGET_NAME}
            PUBLIC $<BUILD_INTERFACE:${BISON_${source_file}_OUTPUT_INCLUDE_DIR}>
        )
    endforeach()
endmacro()
