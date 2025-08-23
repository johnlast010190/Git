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
    Defines the install and package CMake targets.
    Note that install logic for targets is generally added in the
    add_helyx_library() and add_helyx_executable() functions.  If that's not
    possible, the target install logic will be added in the CMakeLists.txt file
    in which the target is defined.

[----------------------------------------------------------------------------]]


# TODO:  Implement debug messaging
message(TITLE "Installation and packing")


# ============================================================================ #
# --------------------------- Binaries Packing ---------------------------- #
# ============================================================================ #

# Install configuration for components binaries
include(packingConfigurationBinaries)


# ============================================================================ #
# ---------------------------- CPack Configuration --------------------------- #
# ============================================================================ #


# setup of variables required for CPackConfig.cmake files generation
if(${HELYX_SYSTEM_NAME} STREQUAL "MSwindows")
    # message(STATUS "Using CPack generator \"ZIP\"")
    set(CPACK_GENERATOR "ZIP")
    set(CPACK_SYSTEM_NAME "win64")
else()
    # message(STATUS "Using CPack generator \"TGZ\"")
    set(CPACK_GENERATOR "TGZ")
    set(CPACK_SYSTEM_NAME "Linux")
endif()

set(CPACK_ARCHIVE_COMPONENT_INSTALL ON)
set(CPACK_PACKAGE_DIRECTORY "${HELYX_PROJECT_DIR}/distribution/")
set(CPACK_CMAKE_GENERATOR "${CMAKE_GENERATOR}")


# ----------------- Binaries and Source Configurations Files ----------------- #

# This is prepended to the component name to create the package file name, e.g.
# ${CPACK_PACKAGE_FILE_NAME}-${COMPONENT_NAME} -> ${CMAKE_PROJECT_NAME}-${HELYX_PROJECT_VERSION}-Windows.zip
# for the core binary package
set(CPACK_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${HELYX_PROJECT_VERSION}")

# message(STATUS "Package name: ${CPACK_PACKAGE_FILE_NAME}-${HELYX_PACKAGE_NAME}")

# Create configuration file for packing binaries
get_cmake_property(CPACK_COMPONENTS_ALL COMPONENTS)
configure_file(
    ${CMAKE_CURRENT_LIST_DIR}/templates/cpackConfig/template-CPackBinariesConfig.cmake
    ${HELYX_PROJECT_DIR}/cbuild/${HELYX_OPTIONS}_${HELYX_MPI_NAME}/CPackBinariesConfig.cmake
    @ONLY
    )

# Create configuration file for packing Core source
include(packingConfigurationCoreSource)

# Create configuration file for packing Modules source
include(packingConfigurationModuleSource)

# Run CPack
# In think we will have to keep include(CPack) to create default config files
# for firt build, even though packing will be called directly from emake
#include(CPack)



# ----------------- Defining commands for packages' targets ----------------- #

# Define custom targets for packing.  This avoids having to deal with package
# targets differently in emake.

# First, remove old-style CPack-generated files.  This is necessary to avoid
# conflicts between the built-in package target and the one defined below.
file(REMOVE
    "${CMAKE_CURRENT_BINARY_DIR}/CPackConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/CPackSourceConfig.cmake"
)

# packageBinaries 'DEPENDS all' only works for CMake version >3.16 and with ninja generator.
# If the requirements are not met, we will directly inform the source directories
if(CMAKE_VERSION VERSION_LESS 3.17)
    message(STATUS "CMake < 3.17 detected, using less robust packageBinaries signature")
    set(packageBinariesDependencies "src-dir;applications-dir;thirdParty-dir")
elseif(CMAKE_GENERATOR MATCHES "Unix Makefiles")
    message(STATUS "Make backend detected, using less robust packageBinaries signature")
    set(packageBinariesDependencies "src-dir;applications-dir;thirdParty-dir")
else()
    set(packageBinariesDependencies "all")
endif()

add_custom_target(package-binaries
    COMMAND ${CMAKE_CPACK_COMMAND} --config ${CMAKE_BINARY_DIR}/CPackBinariesConfig.cmake
    COMMENT "Creating all binary packs, please wait..."
    DEPENDS ${packageBinariesDependencies}
)

add_custom_target(packageBinaries
    COMMAND ${CMAKE_CPACK_COMMAND} --config ${CMAKE_BINARY_DIR}/CPackBinariesConfig.cmake
    COMMENT "Creating all binary packs, please wait..."
    DEPENDS ${packageBinariesDependencies}
)

# For reasons that I don't fully understand, the normal "Make many custom
# targets, then have one that does nothing but depends on all the rest" pattern
# generates the following CPack error:
# CPack Error: Problem removing toplevel directory: HELYXcore-dev/distribution/_CPack_Packages/Linux/TGZ
# I've tried a whole bunch of different ways of doing this, but this is the only
# one that I've found that works.
add_custom_target(packageSource
    DEPENDS dummy_output
    COMMENT "Creating all source packs, please wait..."
)

add_custom_target(package-source
    DEPENDS dummy_output
    COMMENT "Creating all source packs, please wait..."
)

add_custom_command(
    OUTPUT dummy_output
    COMMAND ${CMAKE_CPACK_COMMAND} --config ${CMAKE_BINARY_DIR}/CPackSourcesConfig-HELYX-core.cmake
    COMMENT "Creating HELYXcore source pack, please wait..."
)

# Append commands to the above command (works because they have the same OUTPUT)
set(non_core_modules "${FOUND_MODULES}")
list(REMOVE_ITEM non_core_modules swak4Foam waves2Foam runTimePostProcessing)
foreach(module ${non_core_modules})
    add_custom_command(
        OUTPUT dummy_output APPEND
        COMMAND ${CMAKE_CPACK_COMMAND} --config ${CMAKE_BINARY_DIR}/CPackSourcesConfig-${module}.cmake
        COMMENT "Creating ${module} source packs, please wait..."
    )
endforeach()

# Some versions of CMake (e.g. 3.9, 3.10) reserve the package target name even
# if you dont include(CPack).  Newer versions of CMake (e.g. 3.21) don't do
# this.  We want to always allow setting the package target, but otherwise leave
# the handling of reserved target names untouched.
set(CMAKE_WARN_DEPRECATED FALSE CACHE INTERNAL "")
cmake_policy(SET CMP0037 OLD)
add_custom_target(package
    DEPENDS dummy_output
    COMMAND ${CMAKE_CPACK_COMMAND} --config ${CMAKE_BINARY_DIR}/CPackBinariesConfig.cmake
    COMMENT "Creating all binary packs, please wait..."
    # Modern versions of CMake seem to be able to depend on the "all" target,
    # but we need to do this for backwards-compatibility:
    DEPENDS src-dir applications-dir thirdParty-dir modules-dir
)
cmake_policy(SET CMP0037 NEW)
set(CMAKE_WARN_DEPRECATED TRUE CACHE INTERNAL "")



# ============================================================================ #
# --------------------------------- Exports ---------------------------------- #
# ============================================================================ #

include(CMakePackageConfigHelpers)

set(available_components "${FOUND_MODULES}")
list(REMOVE_ITEM available_components swak4Foam waves2Foam runTimePostProcessing)

# Export core
# message(STATUS "Exporting HELYX-core targets")
install(EXPORT HELYX-core
    DESTINATION ${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}/cmake
    FILE HELYX-core-targets.cmake
    COMPONENT ${HELYX_PACKAGE_NAME}
    )
# Also export targets now for use without installing
export(EXPORT HELYX-core
    FILE ${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}/cmake/HELYX-core-targets.cmake
    )

# Export modules
foreach(comp ${available_components})
    # message(STATUS "Exporting ${comp} targets")
    install(EXPORT ${comp}
        DESTINATION ${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}/cmake
        FILE ${comp}-targets.cmake
        COMPONENT ${HELYX_PACKAGE_NAME}-${comp}
        )
    # Also export targets now for use without installing
    export(EXPORT ${comp}
        FILE ${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}/cmake/${comp}-targets.cmake
        )
endforeach()

install(
    DIRECTORY ${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}/cmake/
    COMPONENT ${HELYX_PACKAGE_NAME}
    DESTINATION ${HELYX_DIR_NAME}/platforms/${HELYX_OPTIONS}/cmake
    PATTERN "*"
    PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE
)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/etc/cmake/templates/template-HELYXConfig.cmake
    "${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}/HELYXConfig.cmake"
    INSTALL_DESTINATION platforms/${HELYX_OPTIONS}
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

# The COMPATIBILITY mode AnyNewerVersion means that the installed package
# version will be considered compatible if it is newer or exactly the same as
# the requested version
write_basic_package_version_file(
    "${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}/HELYXConfigVersion.cmake"
    VERSION 4.4.0
    COMPATIBILITY AnyNewerVersion
)

install(
    FILES
        ${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}/HELYXConfig.cmake
        ${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}/HELYXConfigVersion.cmake
    COMPONENT ${HELYX_PACKAGE_NAME}
    DESTINATION ${HELYX_DIR_NAME}/platforms/${HELYX_OPTIONS}
    PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE
)

message(CLEAN
"[==============================================================================]
")
