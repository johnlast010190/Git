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
    (c) 2021-2022 Engys Ltd.

Description
    Custom findModule for EVTK.

[----------------------------------------------------------------------------]]


# No need to do this, because VTK uses the CONFIG version of find_package
# Push CMAKE_MODULE_PATH to enable calling CMake's find_package() from here
# set(STORED_CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}")
# set(CMAKE_MODULE_PATH "")

set(USING_ENGYS_VTK FALSE)  # Set to true in EVTK config file

# Probably we could get rid of this VTK_REQUIRED check here...
if("${VTK_REQUIRED}" MATCHES "^ON$|^AUTO$")
    if("${VTK_REQUIRED}" STREQUAL "ON")
        set(VTK_REQUIRED_FLAG REQUIRED)
    else()
        set(VTK_REQUIRED_FLAG "")
    endif()

    file(GLOB standard_EVTK_path
        LIST_DIRECTORIES true
        $ENV{HOME}/Engys/HELYX/v*/GUI/ext/
        $ENV{HOME}/StreamlineSolutions/ELEMENTS/v*/GUI/ext/
        )
    list(SORT standard_EVTK_path)
    list(REVERSE standard_EVTK_path)  # Prefer more recent GUI versions

    # For some reason, if VTK_FOUND is set in the cache then the VTK findModule
    # doesn't overwrite it.  Weird.
    unset(VTK_FOUND CACHE)

    find_package(VTK NO_MODULE ${VTK_REQUIRED_FLAG} QUIET
        HINTS
            "${VTK_ARCH_PATH}"
            ${standard_EVTK_path}
            "/opt/ENGYS_VTK/OSMesa_BUILD/"
    )

    set(VTK_FOUND "${VTK_FOUND}" CACHE INTERNAL "")

    # First evaluate the found VTK version to warn and to allow to set
    # VTK_FOUND to FALSE for incompatible versions.
    #
    # For HELYX v4+ we only use VTK 9.2.5, so I'm going to assume that the changes
    # come about because of changes in the major version number rather than
    # anything else.
    # This code could be more elegant...
    if(${VTK_FOUND})
        if("${VTK_VERSION}" VERSION_EQUAL 9.2.5)
            set(VTK_VERSION_9 TRUE CACHE INTERNAL "" FORCE)
        else()
            string(CONCAT s
                "Unrecognised EVTK version \"${VTK_VERSION}\" found at:\n"
                "\t${VTK_DIR}.\n"
                "Officially supported EVTK versions are as follows:\n"
                "\t9.2.5\n"
            )
            if("${VTK_VERSION}" VERSION_LESS 9)
                string(CONCAT s
                    "${s}"
                    "ENGYS VTK was upgraded from version 8.2.0 to 9.2.5 for HELYX "
                    "v4. The found EVTK version is no longer compatible with "
                    "HELYX v4. Therefore, the Run-Time-Post-processing (RTPP) "
                    "will not be compiled. \n"
                    "If RTPP is required, please get the supported EVTK version "
                    "(contact support@engys.com). \n"
                )
                set(VTK_VERSION_9 FALSE CACHE INTERNAL "" FORCE)
            else()
                # probably we will never reach here
                string(CONCAT s
                    "${s}"
                    "HELYX will assume this version is compatible with EVTK 9.2.5, "
                    "and attempt to compile against this EVTK.\n"
                )
                set(VTK_VERSION_9 TRUE CACHE INTERNAL "" FORCE)
            endif()
            message(WARNING "${s}")
        endif()

    endif()

    # Only EVTK version 9 is compatible, the standard VTK is also non compatible
    if(${VTK_VERSION_9})
        if(${USING_ENGYS_VTK})
            string(CONCAT VTK_FOUND_MESSAGE
                "Found EVTK ${VTK_VERSION}:\n"
                "    ENGYS VTK (${VTK_RENDERING_BACKEND}) at ${VTK_DIR}"
            )
        else ()
            set(VTK_FOUND FALSE)
            set(THIRDPARTY_VTK "")
            set(THIRDPARTY_VTK_INC "")
            string(CONCAT VTK_FOUND_MESSAGE
                "Found EVTK ${VTK_VERSION}:\n"
                "    Detected unsupported VTK [not supplied by ENGYS] at: \n"
                "        ${VTK_DIR}.\n"
                "    HELYX supports only the Engys VTK"
            )
        endif()
    else()
        set(VTK_FOUND FALSE)
        set(THIRDPARTY_VTK "")
        set(THIRDPARTY_VTK_INC "")
        string(CONCAT VTK_FOUND_MESSAGE
            "Could NOT find VTK (missing: THIRDPARTY_VTK THIRDPARTY_VTK_INC)\n"
            "    Missing 'VTKConfig' file of a supported EVTK version"
        )
    endif()

else()
    set(VTK_FOUND FALSE)
    string(CONCAT VTK_FOUND_MESSAGE
        "Could NOT find VTK (missing: THIRDPARTY_VTK THIRDPARTY_VTK_INC)\n"
        "    VTK_REQUIRED set to OFF"
    )
endif()

# Expose variables to parent scope
#set(VTK_FOUND "${VTK_FOUND}" CACHE INTERNAL "")
set(VTK_FOUND "${VTK_FOUND}" PARENT_SCOPE)
set(THIRDPARTY_VTK "${THIRDPARTY_VTK}" PARENT_SCOPE)
set(THIRDPARTY_VTK_INC "${THIRDPARTY_VTK_INC}" PARENT_SCOPE)

# This is basically used to print default found/not-found messages
# but (don't know why) its printing all THIRDPARTY_VTK and THIRDPARTY_VTK_INC if VTK is found
#
#   # Handle the QUIETLY and REQUIRED arguments and set VTK_FOUND to TRUE if
#   # all listed variables are TRUE
#   include(FindPackageHandleStandardArgs)
#   find_package_handle_standard_args(VTK DEFAULT_MSG THIRDPARTY_VTK THIRDPARTY_VTK_INC)
message(STATUS "${VTK_FOUND_MESSAGE}")

# No need to do this, because VTK doesn't have a find module!
# Pop CMAKE_MODULE_PATH
# set(CMAKE_MODULE_PATH "${STORED_CMAKE_MODULE_PATH}")
