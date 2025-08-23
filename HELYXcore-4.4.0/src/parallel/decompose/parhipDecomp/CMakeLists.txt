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

[----------------------------------------------------------------------------]]

set(TARGET_NAME "parhipDecomp")

# Don't use add_helyx_library here, because this and OpenFOAM have non-standard
# output directories.  Because it's only two libraries, it's probably not worth
# developing a general method.
add_library(${TARGET_NAME} SHARED
    ./parhipDecomp.C
)

target_include_directories(${TARGET_NAME}
    PRIVATE .
)

add_helyx_include_target(${CMAKE_CURRENT_SOURCE_DIR})
add_to_dir_target(${TARGET_NAME})

helyx_additional_includes(${TARGET_NAME}
    PRIVATE ${THIRDPARTY_MPI_INC}
)

helyx_link_libraries(${TARGET_NAME}
    PUBLIC ${THIRDPARTY_MPI}
    PUBLIC ${THIRDPARTY_PARHIP}
    decompositionMethods
)

# Set output directory (similar to what is done to ptsctochDecomp)
# Note that CMake treats .dll targets as runtime targets
if("${HELYX_SYSTEM_NAME}" STREQUAL "MSwindows")
    set_target_properties(${TARGET_NAME}
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY
        "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${HELYX_MPI_NAME}"
    )
else()
    set_target_properties(${TARGET_NAME}
        PROPERTIES LIBRARY_OUTPUT_DIRECTORY
        "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${HELYX_MPI_NAME}"
    )
endif()

# Manually override installation instructions here - need to change output directory
# For Windows, only pack runtime (.exe and .dll, NOT .dll.a)
if("${HELYX_SYSTEM_NAME}" STREQUAL "MSwindows")
    set(RUNTIME_TOGGLE RUNTIME)
endif()

# This is a core component
install(TARGETS ${TARGET_NAME}
    EXPORT HELYX-core
    ${RUNTIME_TOGGLE}
        DESTINATION ${CMAKE_PROJECT_NAME}-${HELYX_PROJECT_VERSION}/platforms/${HELYX_OPTIONS}/lib/${HELYX_MPI_NAME}
        COMPONENT ${TARGET_PLATFORM_NAME}
        PERMISSIONS
            OWNER_READ OWNER_WRITE OWNER_EXECUTE
            GROUP_READ GROUP_EXECUTE
            WORLD_READ WORLD_EXECUTE
)

# Now that OpenFOAM lives in the $HELYX_MPI_NAME subdir, we need to add
# $ORIGIN/.. and $ORIGIN to the RUNPATH in order for run-time selection to work
get_target_property(install_runpath ${TARGET_NAME} INSTALL_RPATH)
string(REPLACE "$ORIGIN/../lib"
    "$ORIGIN;$ORIGIN/.." install_runpath
    "${install_runpath}"
)
set_target_properties(${TARGET_NAME}
    PROPERTIES INSTALL_RPATH "${install_runpath}"
)

# Append to list of Helyx targets for emake autocompletion
get_property(TARGETS GLOBAL PROPERTY ALL_HELYX_LIB_TARGETS)
list(APPEND TARGETS ${TARGET_NAME})
set_property(GLOBAL PROPERTY ALL_HELYX_LIB_TARGETS ${TARGETS})
