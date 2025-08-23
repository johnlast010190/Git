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
    Defines package configuration file of modules.

[----------------------------------------------------------------------------]]


# ============================================================================ #
# -------------------------- Module Sources Packing -------------------------- #
# ============================================================================ #

# Removing swak4Foam, waves2Foam and runTimePostProcessing from HELYX add-ons list
set(available_modules "${FOUND_MODULES}")
list(REMOVE_ITEM available_modules swak4Foam waves2Foam runTimePostProcessing)

# Create source configuration file for HELYX add-ons
foreach(module_name ${available_modules})

    if(NOT "${module_name}" IN_LIST FOUND_MODULES)
        continue()
    endif()

    set(CPACK_IGNORE_FILES "")

    # Note: CPACK_INSTALLED_DIRECTORIES is a list of pairs
    set(CPACK_INSTALLED_DIRECTORIES
        "${HELYX_PROJECT_DIR}/modules/${module_name}"  # Dir to install
        "${HELYX_DIR_NAME}/modules/${module_name}" # Relative location in package
    )

    # message(STATUS "Packing module ${module_name} sources")

    set(module_whitelist
        "src;tutorials;CMakeLists.txt;applications"
    )

    # Each module may have slightly different contents
    if ("${module_name}" STREQUAL "HELYX-adjoint")
        list(APPEND module_whitelist bin)
    elseif ("${module_name}" STREQUAL "HELYX-coupled")
    elseif ("${module_name}" STREQUAL "HELYX-hydro")
    elseif ("${module_name}" STREQUAL "HELYX-marine")
    elseif ("${module_name}" IN_LIST HELYX_ADDITIONAL_MODULES)
        set(module_whitelist "")
    endif()

    whitelist_items_in_dir("${HELYX_PROJECT_DIR}/modules/${module_name}"
        "${module_whitelist}"
    )

    # Creating cpack configuration file for modules source
    set(CPACK_PRE_BUILD_SCRIPTS "")
    set(CPACK_PACKAGE_FILE_NAME
        "${HELYX_DIR_NAME}-${HELYX_PACKAGE_NAME}-${module_name}-SRC"
    )
    configure_file(
        ${CMAKE_CURRENT_LIST_DIR}/templates/cpackConfig/template-CPackSourcesConfig.cmake
        ${HELYX_PROJECT_DIR}/cbuild/${HELYX_OPTIONS}_${HELYX_MPI_NAME}/CPackSourcesConfig-${module_name}.cmake
        @ONLY
    )

    match_installed_source_to_packed_source()

endforeach()



