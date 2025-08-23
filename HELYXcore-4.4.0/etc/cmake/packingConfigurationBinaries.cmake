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
    (c) 2019 Engys Ltd.

Description
    Defines the install and package of binaries.

[----------------------------------------------------------------------------]]


# ============================================================================ #
# --------------------------- Core Binary Packing ---------------------------- #
# ============================================================================ #

# message(STATUS "Packing HELYX-core")

install(
    DIRECTORY bin etc
    COMPONENT ${HELYX_PACKAGE_NAME}
    DESTINATION ${HELYX_DIR_NAME}
    USE_SOURCE_PERMISSIONS
)

# For some reason, the emake permissions appear to be going missing...
# Installing emake explicitly appears to fix the problem...
install(
    FILES bin/emake
    COMPONENT ${HELYX_PACKAGE_NAME}
    DESTINATION ${HELYX_DIR_NAME}/bin
    PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE
)

# User settings file may be a symlink, in which case install a copy
# If userSettings.cmake is a symlink, it keep the symlink and the real file
# is not copied to the package, but the real file is indeed copied over for installation
get_filename_component(real_settings_file ${HELYX_SETTINGS_FILE} REALPATH)
install(
    FILES ${real_settings_file}
    COMPONENT ${HELYX_PACKAGE_NAME}
    DESTINATION ${HELYX_DIR_NAME}/etc
    RENAME userSettings.cmake
    PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE
)

install(
    FILES LICENSE README
    COMPONENT ${HELYX_PACKAGE_NAME}
    DESTINATION ${HELYX_DIR_NAME}
)

install(
    DIRECTORY examples
    COMPONENT ${HELYX_PACKAGE_NAME}
    DESTINATION ${HELYX_DIR_NAME}
    USE_SOURCE_PERMISSIONS
)


# Installing additional files required by Windows, including ThirdParty
if(${HELYX_SYSTEM_NAME} STREQUAL "MSwindows")

    # message(STATUS
    #     "The following files will be written to facilitate HELYXcore in Windows:
    #     HELYX_Terminal.bat
    #     activeBuild.bat
    #     ${HELYX_OPTIONS}.bat"
    # )

    configure_file(
        ${CMAKE_CURRENT_LIST_DIR}/templates/template-batchrc.bat
        ${HELYX_PROJECT_DIR}/platforms/${HELYX_OPTIONS}.bat
        @ONLY
        NEWLINE_STYLE WIN32
    )

    install(
        FILES platforms/${HELYX_OPTIONS}.bat
        COMPONENT ${HELYX_PACKAGE_NAME}
        DESTINATION ${HELYX_DIR_NAME}/platforms/
        PERMISSIONS
            OWNER_READ OWNER_WRITE OWNER_EXECUTE
            GROUP_READ GROUP_EXECUTE
            WORLD_READ WORLD_EXECUTE
    )

    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/HELYX_Terminal.bat
        "@echo off \r\n\r\ncmd /k \"%~dp0\\${HELYX_DIR_NAME}\\platforms\\activeBuild-${HELYX_PRECISION_OPTION}.bat\""
    )

    install(
        FILES ${CMAKE_CURRENT_BINARY_DIR}/HELYX_Terminal.bat
        COMPONENT ${HELYX_PACKAGE_NAME}
        DESTINATION .
    )

    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/activeBuild.bat
        "@echo off \r\n\r\ncall \"%~dp0\\${HELYX_OPTIONS}.bat\""
    )

    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/activeBuild-${HELYX_PRECISION_OPTION}.bat
        "@echo off \r\n\r\ncall \"%~dp0\\${HELYX_OPTIONS}.bat\""
    )

    install(
        FILES
            ${CMAKE_CURRENT_BINARY_DIR}/activeBuild.bat
            ${CMAKE_CURRENT_BINARY_DIR}/activeBuild-${HELYX_PRECISION_OPTION}.bat
        COMPONENT ${HELYX_PACKAGE_NAME}
        DESTINATION ${HELYX_DIR_NAME}/platforms
        PERMISSIONS
            OWNER_READ OWNER_WRITE OWNER_EXECUTE
            GROUP_READ GROUP_EXECUTE
            WORLD_READ WORLD_EXECUTE
        )

    # Pack ThirdParty alongside Core
    install(
        DIRECTORY ${HELYX_THIRDPARTY_DIR}/
        COMPONENT ${HELYX_PACKAGE_NAME}
        DESTINATION ThirdParty-${HELYX_THIRDPARTY_VERSION}
    )

    # Pack wtee and wtail
    # message(STATUS "Packing wtee.exe and wtail.exe")
    install(
        FILES ${HELYX_THIRDPARTY_DIR}/wtee.exe ${HELYX_THIRDPARTY_DIR}/wtail.exe
        COMPONENT ${HELYX_PACKAGE_NAME}
        DESTINATION ${HELYX_DIR_NAME}/platforms/${HELYX_OPTIONS}/bin
        PERMISSIONS
            OWNER_READ OWNER_WRITE OWNER_EXECUTE
            GROUP_READ GROUP_EXECUTE
            WORLD_READ WORLD_EXECUTE
    )
    # May as well do this here as well...
    # This is a horrible hack, but it's what's done in the Allwmake file for Windows.
    file(COPY "${HELYX_THIRDPARTY_DIR}/wtee.exe" DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
    file(COPY "${HELYX_THIRDPARTY_DIR}/wtail.exe" DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

else()

    # Pack linux sourceable files
    install(
        FILES
            platforms/activeBuild.shrc
            platforms/unsetActiveBuild.shrc
            platforms/${HELYX_OPTIONS}_${HELYX_MPI_NAME}.shrc
            platforms/activeBuild-${HELYX_PRECISION_OPTION}.shrc
        COMPONENT ${HELYX_PACKAGE_NAME}
        DESTINATION ${HELYX_DIR_NAME}/platforms
        PERMISSIONS
            OWNER_READ OWNER_WRITE OWNER_EXECUTE
            GROUP_READ GROUP_EXECUTE
            WORLD_READ WORLD_EXECUTE
    )

    # Pack linux extra-sourceable files
    install(
        DIRECTORY cbuild/${HELYX_OPTIONS}_${HELYX_MPI_NAME}/extraSourceableFiles
        COMPONENT ${HELYX_PACKAGE_NAME}
        DESTINATION ${HELYX_DIR_NAME}/platforms/${HELYX_OPTIONS}/cmake
        USE_SOURCE_PERMISSIONS
    )

    # Pack additional mpi-dependen libraries if requested
    if(NOT "" STREQUAL "${HELYX_PACK_EXTRA_MPI_DEPENDEND_LIBS}")
        install(
            FILES
                platforms/${HELYX_OPTIONS}_${HELYX_PACK_EXTRA_MPI_DEPENDEND_LIBS}.shrc
            COMPONENT ${HELYX_PACKAGE_NAME}
            DESTINATION ${HELYX_DIR_NAME}/platforms
            PERMISSIONS
                OWNER_READ OWNER_WRITE OWNER_EXECUTE
                GROUP_READ GROUP_EXECUTE
                WORLD_READ WORLD_EXECUTE
        )
        install(
            DIRECTORY platforms/${HELYX_OPTIONS}/lib/${HELYX_PACK_EXTRA_MPI_DEPENDEND_LIBS}
            COMPONENT ${HELYX_PACKAGE_NAME}
            DESTINATION ${HELYX_DIR_NAME}/platforms/${HELYX_OPTIONS}/lib
            USE_SOURCE_PERMISSIONS
            FILE_PERMISSIONS
                OWNER_READ OWNER_WRITE OWNER_EXECUTE
                GROUP_READ GROUP_EXECUTE
                WORLD_READ WORLD_EXECUTE
        )
    endif()

    # Pack ThirdParty alongside Core for POSIX
    install(
        DIRECTORY ${HELYX_THIRDPARTY_DIR}/platforms
        COMPONENT ${HELYX_PACKAGE_NAME}
        DESTINATION ThirdParty-${HELYX_PROJECT_VERSION}
        PATTERN "*"
        PERMISSIONS
            OWNER_READ OWNER_WRITE OWNER_EXECUTE
            GROUP_READ GROUP_EXECUTE
            WORLD_READ WORLD_EXECUTE
    )

    # Write .info at install time and pack alongside core
    install(CODE "file(WRITE ${HELYX_PROJECT_DIR}/.info
    \"${HELYX_PACKAGE_NAME}\"
    )")
    install(FILES .info COMPONENT ${HELYX_PACKAGE_NAME} DESTINATION .)

    install(
        FILES ${HELYX_THIRDPARTY_DIR}/DISCLAIMER
        COMPONENT ${HELYX_PACKAGE_NAME}
        DESTINATION ThirdParty-${HELYX_PROJECT_VERSION}
    )

    # We no longer build libOceanWave3D
    # libOceanWave3D is a special case - it's built from a makefile in
    # waves2Foam/ThirdParty.  Install logic therefore can't be imparted at build
    # time, so we have to treat it here.
    # if(waves2Foam IN_LIST FOUND_MODULES)
    #     install(
    #         FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libOceanWave3D.so
    #         COMPONENT ${HELYX_PACKAGE_NAME}
    #         DESTINATION ${HELYX_DIR_NAME}/platforms/${HELYX_OPTIONS}/lib
    #         PERMISSIONS
    #             OWNER_READ OWNER_WRITE OWNER_EXECUTE
    #             GROUP_READ GROUP_EXECUTE
    #             WORLD_READ WORLD_EXECUTE
    #     )
    # endif()

endif()


# ============================================================================ #
# --------------------------- Module Binary Packing -------------------------- #
# ============================================================================ #

# ---------------------------------  Macros  --------------------------------- #

# Module name should be something like HELYX-adjoint
# Don't use this for Adjoint plus
macro(pack_module_bin module_name)
    install(
        DIRECTORY modules/${module_name}/tutorials
        COMPONENT ${HELYX_PACKAGE_NAME}-${module_name}
        DESTINATION ${HELYX_DIR_NAME}/modules/${module_name}/
        PATTERN "*"
        PERMISSIONS
            OWNER_READ OWNER_WRITE OWNER_EXECUTE
            GROUP_READ GROUP_EXECUTE
            WORLD_READ WORLD_EXECUTE
    )
endmacro()


# --------------------------  Other module packing  -------------------------- #

foreach(module IN ITEMS HELYX-adjoint HELYX-coupled HELYX-marine HELYX-hydro ${HELYX_ADDITIONAL_MODULES})
    if("${module}" IN_LIST FOUND_MODULES)
        # message(STATUS "Packing module ${module}")
        pack_module_bin(${module})
    endif()
endforeach()
