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

# This is used occasionally in CMakeLists files, may as well leave it here
set(WAVES_SRC "${CMAKE_CURRENT_LIST_DIR}/src")

# GSL removed in version 3.4.
# This probably isn't the most elegant way of doing things - this function only
# works like this for historical reasons. Consider refactoring.
thirdParty_package_deprecation_warning("GSL" "3.4")
if (NOT "" STREQUAL "${GSL_FOUND_MESSAGE}")
    if ("ON" STREQUAL GSL_REQUIRED)
        set(MISSING_REQUIRED_LIB TRUE)
        message(SEND_ERROR "${GSL_FOUND_MESSAGE}")
    else()
        message(STATUS "${GSL_FOUND_MESSAGE}")
    endif()
endif()

add_definitions(-DOFPLUSBRANCH=1)

# Note that this evaluates to the OpenFOAM version HELYX is based on (e.g. 1706)
# It should NOT evaluate to the HELYX vesion (e.g. dev, 3.2.0, etc...)
add_definitions(-DOFVERSION=1706)

add_subdirectory(applications)

add_subdirectory(src)
