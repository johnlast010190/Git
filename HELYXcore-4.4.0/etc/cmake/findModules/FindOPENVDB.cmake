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
    (c) 2022 Engys Ltd.

Description
    Custom findModule for OpenVDB

[----------------------------------------------------------------------------]]


# N.B.  OpenVDB headers are not in standard include dir, instead, they are located
# in */include/openvdb/, but OpenVDB own code expects includes from */include/.
# Inside 'helyx_find_thirdparty_package' the found includes path is corrected to
# point to */include/ (instead of */include/openvdb/) and then, the headers are
# associated to the OpenVDB libs.
# This is done by listing the headers in 'helyx_find_thirdparty_package' as a
# suffix of the expect include folder.
helyx_find_thirdparty_package(OPENVDB "openvdb" "openvdb/openvdb.h")

# If OPENVDB compilation is disabled, simply return from this script
# There is no need to evaluate dependencies
if(OPENVDB_REQUIRED STREQUAL "OFF" OR NOT OPENVDB_FOUND)
    return()
endif()

## N.B.  Could add extra checks and detection here to not have BOOST and CBLOSC
## on OPTIONAL_THIRDPARTY_LIBRARIES:
## Check combinations of OPENVDB_REQUIRED and dependency_REQUIRED
## then set ${OPENVDB|dependency}_FOUND or ${OPENVDB|dependency}_REQUIRED accordingly
#helyx_check_module_mandatory_dependency(OPENVDB TBB)
#helyx_check_module_mandatory_dependency(OPENVDB BOOST)
#helyx_check_module_mandatory_dependency(OPENVDB CBLOSC)

# OpenVDB dependes on TBB, Boost and CBLOSC.
# Link the dependencies if they are found, otherwise, disable OpenVDB and send a
# warning message. At this point, all dependencies should have been evaluated.
string(CONCAT warn_message_libs
"OPENVDB libraries found, but dependencies are missing:
    TBB_FOUND: ${TBB_FOUND}
    BOOST_FOUND: ${BOOST_FOUND}
    CBLOSC_FOUND: ${CBLOSC_FOUND}
Review the dependencies setup in the settings file and reconfigure.
To enable <depenency>, check if:
    <depenency> is built in ThirdParty or available in the system,
    <depenency>_ARCH_PATH is properly set in the userSettings file,
    <depenency>_REQUIRED is set to ON|AUTO.\n
- Setting OPENVDB_REQUIRED to OFF.
- Setting OPENVDB_FOUND to FALSE.\n"
)
if (${OPENVDB_FOUND})
    if(${TBB_FOUND})
        target_link_libraries("openvdb" INTERFACE ${THIRDPARTY_TBB})
    endif()
    if(${BOOST_FOUND})
        target_link_libraries("openvdb" INTERFACE ${THIRDPARTY_BOOST})
    endif()
    if(${CBLOSC_FOUND})
        target_link_libraries("openvdb" INTERFACE ${THIRDPARTY_CBLOSC})
    endif()
    if ("FALSE" STREQUAL "${TBB_FOUND}" OR
        "FALSE" STREQUAL "${BOOST_FOUND}" OR
        "FALSE" STREQUAL "${CBLOSC_FOUND}")
        # Disable OpenVDB
        set(OPENVDB_REQUIRED OFF CACHE INTERNAL "")
        set(OPENVDB_FOUND FALSE CACHE INTERNAL "")
        # Overwrite default message
        set("OPENVDB_FOUND_MESSAGE" "OPENVDB  (OPENVDB_REQUIRED set to OFF due to missing dependencies)" )
        # print warning
        message(WARNING "${warn_message_libs}")
    endif()
endif()


# BUG:  If OpenVDB was compiled agains Boost and C-Blosc, but their *_ARCH_PATH
# are not set in userSettings.cmake there will be a linking issue on helyxVDB
# In this situation, Boost and C-Blosc wont be found, and their libraries will
# be linked to their absolute paths. If the binaries are copied over another machine,
# these absolute paths wont exist...
#
# TO DO:  Properly deal with thirdParty_<package> dependencies:
# need to create a funtion to standardize that, which should also check
# all possible combinations of <dependency>_REQUIRED and <dependency>_ARCH_PATH
#
# ps.:  helyx_find_thirdparty_package needs that <package>_REQUIRED is defined, otherwise
# it rise an error. We could add BOOST and CBLOSC to OPTIONAL_THIRDPARTY_LIBRARIES,
# but then, a Find<package> will need to be created. Since, BOOST and CBLOSC are
# dependencies of OpenVDB, we should keed their finding here.
