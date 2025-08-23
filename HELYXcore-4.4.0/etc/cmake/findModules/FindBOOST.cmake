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
    Custom findModule for TBB

[----------------------------------------------------------------------------]]


# N.B.  Boost headers are not in standard include dir, instead, they are located
# in */include/boost/, but Boost own code expects includes from */include/.
# Inside 'helyx_find_thirdparty_package' the found includes path is corrected to
# point to */include/ (instead of */include/boost/) and then, the headers are
# associated to the BOOST libs.
# This is done by listing the headers in 'helyx_find_thirdparty_package' as a
# suffix of the expect include folder.
helyx_find_thirdparty_package(BOOST "boost_iostreams" "boost/config.hpp")


## I THINK WE CAN REMOVE THE CODE BELOW:
#
## CMake may still find system boost outside the BOOST_ARCH_PATH,
## e.g. if BOOST_ARCH_PATH is defined but the path doesn`t exist.
## The found boost path and the BOOST_ARCH_PATH need to be in accordance.
## Otherwise, set BOOST_FOUND FALSE.
## Also, things like SloanRenumber uses BOOST_ARCH_PATH for linking and includes.
#string(CONCAT warn_message_libs
#"BOOST libraries found, but not on BOOST_ARCH_PATH
#    BOOST_LIBRARIES:
#        \"${THIRDPARTY_BOOST_DIR}\"
#    BOOST_INCLUDES:
#        \"${THIRDPARTY_BOOST_INC}\"
#    BOOST_ARCH_PATH:
#        \"${BOOST_ARCH_PATH}\"
#If these are the correct libraries, consider changing BOOST_ARCH_PATH.
#Setting BOOST_FOUND to FALSE\n"
#)
#
#if(BOOST_FOUND)
#    compare_to_arch_path("BOOST" "${THIRDPARTY_BOOST_DIR}")
#    if (NOT MATCH_BOOST_ARCH_PATH)
#        set(BOOST_FOUND FALSE CACHE INTERNAL "")
#        set(BOOST_REQUIRED OFF CACHE INTERNAL "")
#        message(WARNING "${warn_message_libs}")
#        set("BOOST_FOUND_MESSAGE" "BOOST  (BOOST_REQUIRED set to OFF due to BOOST_ARCH_PATH missmatch)" )
#    endif()
#endif()
