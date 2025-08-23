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

Description
    Defines package configuration file of sources.

[----------------------------------------------------------------------------]]



function(whitelist_items_in_dir dir whitelist_items)
    if(NOT ${ARGC} EQUAL 2)
        string(CONCAT s
            "whitelist_items_in_dir() requires exactly two arguments, "
            "but ${ARGC} were provided"
        )
        message(FATAL_ERROR "${s}")
    endif()

    file(GLOB blacklist RELATIVE "${dir}" "${dir}/*")

    # There's a nasty edge-case where a file and a directory have the same name
    # Should really deal with this properly, but for now we'll just warn
    set(test "${blacklist}")
    list(REMOVE_DUPLICATES test)
    if(NOT "${test}" STREQUAL "${blacklist}")
        string(CONCAT s
            "Got duplicates in blacklist for \"${dir}\":\n"
            "\t${blacklist}"
            "Your source package may not be what you're expecting"
        )
        message(WARNING "${s}")
    endif()

    foreach(whitelist_item ${whitelist_items})
        if ("${whitelist_item}" IN_LIST blacklist)
            list(REMOVE_ITEM blacklist "${whitelist_item}")
        else()
            string(CONCAT s
                "Could not find \"${whitelist_item}\" in \"${dir}\" in order "
                "to whitelist it for packing"
            )
            message(AUTHOR_WARNING "${s}")
        endif()
    endforeach()

    get_filename_component(dir_name "${dir}" NAME)
    foreach(blacklist_item ${blacklist})
        if(IS_DIRECTORY ${dir}/${blacklist_item})
            # message("Blacklisting directory \"${dir_name}/${blacklist_item}/\"")
            # Must append / to dir names, to avoid accidentally catching files
            list(APPEND CPACK_IGNORE_FILES "${dir_name}/${blacklist_item}/")
        else()
            # Can't check for existance here - might be a broken symlink!
            list(APPEND CPACK_IGNORE_FILES "${dir_name}/${blacklist_item}")
        endif()
    endforeach()

    set(CPACK_IGNORE_FILES "${CPACK_IGNORE_FILES}" PARENT_SCOPE)
endfunction()


# This macro uses CPACK_IGNORE_FILES and CPACK_INSTALLED_DIRECTORIES to call
# install() in a way that matches the source packs that CPack generates.
macro(match_installed_source_to_packed_source)
    set(install_ignore_files "${CPACK_IGNORE_FILES}")

    # It's unnecessarily awkward to adds exclude patterns to install()
    # programatically, so we'll just make it one massive regex...
    string(REGEX REPLACE "/?;" "$|" ignore_string "(${install_ignore_files})")

    # Need to install sources in order for external projects to have access to
    # code!
    list(LENGTH CPACK_INSTALLED_DIRECTORIES len)
    math(EXPR len "${len}-1")
    foreach(i RANGE 0 ${len} 2)
        list(GET CPACK_INSTALLED_DIRECTORIES ${i} dir)
        math(EXPR i+1 "${i}+1")
        list(GET CPACK_INSTALLED_DIRECTORIES ${i+1} destination)
        install(DIRECTORY "${dir}/"
            DESTINATION "./${destination}"
            USE_SOURCE_PERMISSIONS
            COMPONENT SOURCE
            REGEX "${ignore_string}" EXCLUDE
        )
    endforeach()
endmacro()



# ============================================================================ #
# ------------------------ Files for Core Source Pack ------------------------ #
# ============================================================================ #

# message(STATUS "Packing ${CMAKE_PROJECT_NAME} sources")

# Note: CPACK_INSTALLED_DIRECTORIES is a list of pairs:
#   ${dir_to_install};${relative_location_in_package}

# Blacklisting unnecessary files/folders for ThirdParty
if(${HELYX_SYSTEM_NAME} STREQUAL "POSIX")
    get_filename_component(THIRDPARTY_DIR_NAME "${HELYX_THIRDPARTY_DIR}" NAME)
    list(APPEND CPACK_INSTALLED_DIRECTORIES "${HELYX_THIRDPARTY_DIR};${THIRDPARTY_DIR_NAME}")

    set(thirdParty_whitelist "CMakeLists.txt;DISCLAIMER;etc")

    if (EXISTS ${HELYX_THIRDPARTY_DIR}/downloads)
        message(STATUS "Packing third party downloads directory")
        list(APPEND thirdParty_whitelist "downloads")
    else()
        message(STATUS
            "Not packing third party archives (\\\${HELYX_THIRDPARTY_DIR}/downloads does not exist)"
            )
    endif()

    if (EXISTS ${HELYX_THIRDPARTY_DIR})
        whitelist_items_in_dir("${HELYX_THIRDPARTY_DIR}" "${thirdParty_whitelist}")
    else()
        message(STATUS
            "Not packing third party (HELYX_THIRDPARTY_DIR does not exist)"
            )
    endif()

endif()

list(APPEND CPACK_INSTALLED_DIRECTORIES "${HELYX_PROJECT_DIR};${HELYX_DIR_NAME}")
list(APPEND core_whitelist
    "CMakeLists.txt"
    "LICENSE"
    "README"
    "applications"
    "bin"
    "doc"
    "etc"
    "modules"
    "src"
    "thirdParty"
    "examples"
    "VTK"
)
whitelist_items_in_dir("${HELYX_PROJECT_DIR}" "${core_whitelist}")

list(APPEND modules_whitelist
    "CMakeLists.txt"
    )
foreach(core_module waves2Foam;swak4Foam;runTimePostProcessing)
    if("${core_module}" IN_LIST FOUND_MODULES)
        list(APPEND modules_whitelist "${core_module}")
    endif()
endforeach()
whitelist_items_in_dir("${HELYX_PROJECT_DIR}/modules" "${modules_whitelist}")

# CMake insists in packing "HELYXcore-dev/distribution/_CPack_Packages/Linux/TGZ"
# We can't include /distribution/ because there are other folders called distribution inside src
# Therefore, we remove _CPack_Packages instead, leving the cpack distribution folder empty
#
# p.s.:  It looks like CMake doesn`t like to remove *.cmake files. It accepts 'userSettings.txt'
# but not 'userSettings.cmake'
list(APPEND CPACK_IGNORE_FILES
    "/_CPack_Packages/"
    ${CMAKE_INSTALL_PREFIX}
)

match_installed_source_to_packed_source()

# N.B.  We can`t use install(FILE ...) to copy the etc/exampleLinuxSettings.cmake file
# to etc/userSettings.cmake for the package because it only do the copy for the
# installation step and not for the packing step.



# ============================================================================ #
# --------------------------- Generate CPack Files --------------------------- #
# ============================================================================ #

# Creating cpack configuration file for core source
# "module_name" and "CPACK_PRE_BUILD_SCRIPTS" are required in template-CPackSourcesConfig.cmake
set(CPACK_PRE_BUILD_SCRIPTS "${HELYX_PROJECT_DIR}/cbuild/${HELYX_OPTIONS}_${HELYX_MPI_NAME}/CPackSymlinksSource.cmake")
set(CPACK_PACKAGE_FILE_NAME "${HELYX_DIR_NAME}-${HELYX_PACKAGE_NAME}-SRC")
configure_file(
    ${CMAKE_CURRENT_LIST_DIR}/templates/cpackConfig/template-CPackSourcesConfig.cmake
    ${HELYX_PROJECT_DIR}/cbuild/${HELYX_OPTIONS}_${HELYX_MPI_NAME}/CPackSourcesConfig-HELYX-core.cmake
    @ONLY
    )

# Creating cpack configuration file for symlinks in core
# "real_settings_file" and "cpack_temporary_source_dir" are required in templates-installSymlinks.cmake
if(EXISTS "${HELYX_PROJECT_DIR}/examples")

    get_filename_component(real_settings_file ${HELYX_SETTINGS_FILE} REALPATH)

    set(cpack_temporary_source_dir
        ${CPACK_PACKAGE_DIRECTORY}_CPack_Packages/${CPACK_SYSTEM_NAME}/${CPACK_GENERATOR}/${HELYX_DIR_NAME}-${HELYX_PACKAGE_NAME}-SRC/${HELYX_DIR_NAME}/
        )

    configure_file(
        ${CMAKE_CURRENT_LIST_DIR}/templates/cpackConfig/templates-installSymlinks.cmake
        ${HELYX_PROJECT_DIR}/cbuild/${HELYX_OPTIONS}_${HELYX_MPI_NAME}/CPackSymlinksSource.cmake
        @ONLY
        )
endif()
