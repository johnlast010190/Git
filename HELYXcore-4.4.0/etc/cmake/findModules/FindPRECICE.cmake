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
    (c) 2023 Engys Ltd.

Description
Custom findModule for preCICE

[----------------------------------------------------------------------------]]


# N.B.  PreCICE headers are not in standard include dir, instead, they are located
# in */include/precise/, but PreCICE own code expects includes from */include/.
# Inside 'helyx_find_thirdparty_package' the found includes path is corrected to
# point to */include/ (instead of */include/precise/) and then, the headers are
# associated to the PreCICE libs.
# This is done by listing the headers in 'helyx_find_thirdparty_package' as a
# suffix of the expect include folder.
helyx_find_thirdparty_package(PRECICE
    "precice"
    "precice/SolverInterface.hpp"
)

# If PRECICE compilation is disabled, simply return from this script.
# There is no need to evaluate dependencies
if(PRECICE_REQUIRED STREQUAL "OFF" OR NOT PRECICE_FOUND)
    return()
endif()

# N.B.  PreCICE requires boost, Eigen3 and LibXML2
# Eigen3 and LibXML2 should be installed in the system, so does it make sense
# to add find_package() for two libs?

## N.B.  Could add extra checks and detection here to not have BOOST on
## OPTIONAL_THIRDPARTY_LIBRARIES:
## Check combinations of PRECISE_REQUIRED and dependency_REQUIRED
## then set ${PRECISE|dependency}_FOUND or ${PRECISE|dependency}_REQUIRED accordingly
#helyx_check_module_mandatory_dependency(PRECISE BOOST)

# Precise depends on BOOST.
# Link BOOST if found, otherwise, disable PRECISE and send a warning message.
# At this point, all dependencies should have been evaluated.
STRING(CONCAT warn_message_libs
"PRECISE libraries found, but optional dependencies are missing:
    BOOST_FOUND: ${BOOST_FOUND}
If you are willing to use BOOST, review the dependencies setup in the settings file and reconfigure.
To enable 'BOOST', check if:
    BOOST is built in ThirdParty or available in the system,
    BOOST_ARCH_PATH is properly set in the userSettings file,
    BOOST_REQUIRED is set to ON|AUTO.\n"
)
if (${PRECISE_FOUND})
    if(${BOOST_FOUND})
        target_link_libraries("precice" INTERFACE ${THIRDPARTY_BOOST})
    else()
        # Disable OpenCascade
        set(PRECISE_REQUIRED OFF CACHE INTERNAL "")
        set(PRECISE_FOUND FALSE CACHE INTERNAL "")
        # Ovewrite default message
        set("PRECISE_FOUND_MESSAGE" "PRECISE  (PRECISE_REQUIRED set to OFF due to missing dependencies)" )
        message(WARNING "${warn_message_libs}")
    endif()
endif()
