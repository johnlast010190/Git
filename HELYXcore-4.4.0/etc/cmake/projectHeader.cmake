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
    Include file for the top-level CMakeLists.txt in a project. Can be used
    in a custom project outside the source tree.

[----------------------------------------------------------------------------]]


# N.B.: EMAKE_BUILD=TRUE when we ran "emake -r", so it won`t run the script
# below. But at the at projectFooter.cmake, we set EMAKE_BUILD=FALSE (to allow
# configuration change from tools like VSCode). Now, EMAKE_BUILD is FALSE when
# "emake" runs, triggering the script below and re-loading the userSettings.
# If inner CMake variables (like CMAKE_*_COMPILER) are set in userSettings, they
# will overwritte the values computed by the automatic system-detection,
# changing the build rules and forcing a full rebuild.
#
# To avoid this, EMAKE_BUILD=FALSE could be set in userSettings and removed from
# the projectFooter.cmake, but that may cause too much confusion...
# Best way for now is simply ensure that no CMake variables are overwritten
# by the userSettings. So far only CMAKE_*_COMPILER matters.
#
# Note: CMake does not allow changing compilers without completely regenerating
# the buildsystem, so the above appears to be an effort to fix something that is
# not fixable (or particularly desirable to fix: just use multiple object
# directories!)
if(NOT EMAKE_BUILD)
    # Load userSettings and generate minimum sourceable files
    # message(STATUS "Not an emake build, calling emake pre-requisites as a script")
    include("${CMAKE_CURRENT_LIST_DIR}/scripts/generateEmakePrerequisites.cmake")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/functions.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/userSettingsConfiguration.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/cmakeConfiguration.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/thirdPartyConfiguration.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/buildConfiguration.cmake")
