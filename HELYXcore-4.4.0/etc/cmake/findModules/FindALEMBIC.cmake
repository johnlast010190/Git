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
    Custom findModule for ALEMBIC

[----------------------------------------------------------------------------]]


# N.B.  Alembic headers are not in standard include dir, instead, they are located
# in */include/Alembic/**/, but Alembic own code expects includes from */include/
# Inside 'helyx_find_thirdparty_package' the found includes path is corrected to
# point to */include/ (instead of */include/Alembic/) and then, the headers are
# associated to the Alembic libs.
# This is done by listing the headers in 'helyx_find_thirdparty_package' as a
# suffix of the expect include folder.
set(alembic_headers
    "Alembic/Abc/All.h"
    "Alembic/AbcCollection/All.h"
    "Alembic/AbcCoreAbstract/All.h"
    "Alembic/AbcCoreFactory/All.h"
    "Alembic/AbcCoreLayer/Util.h"
    "Alembic/AbcCoreOgawa/All.h"
    "Alembic/AbcGeom/All.h"
    "Alembic/AbcMaterial/All.h"
    "Alembic/Util/All.h"
    )

## 'helyx_find_thirdparty_package' looks for includes in lower-case suffixes.
## Because of uppercase 'A' of Alembic, we need to set ALEMBIC_PATH_SUFFIXES
#set(ALEMBIC_PATH_SUFFIXES "include/Alembic")
helyx_find_thirdparty_package(ALEMBIC "Alembic" "${alembic_headers}")

# If ALEMBIC compilation is disabled, simply return from this script.
# There is no need to evaluate dependencies
if(ALEMBIC_REQUIRED STREQUAL "OFF" OR NOT ALEMBIC_FOUND)
    return()
endif()

## N.B.  Could add extra checks and detection here to not have IMATH on
## OPTIONAL_THIRDPARTY_LIBRARIES:
## Check combinations of ALEMBIC_REQUIRED and dependency_REQUIRED
## then set ${ALEMBIC|dependency}_FOUND or ${ALEMBIC|dependency}_REQUIRED accordingly
#helyx_check_module_mandatory_dependency(ALEMBIC IMATH)

# Link IMath if found, otherwise, disable Alembic and send a warning message.
# At this point, all dependencies should have been evaluated.
string(CONCAT warn_message_libs
"ALEMBIC libraries found, but dependencies are missing:
    IMATH_FOUND: ${IMATH_FOUND}
Review the dependencies setup in the settings file and reconfigure.
To enable 'IMATH', check if:
    IMATH is built in ThirdParty or available in the system,
    IMATH_ARCH_PATH is properly set in the userSettings file,
    IMATH_REQUIRED is set to ON|AUTO.\n
- Setting ALEMBIC_REQUIRED to OFF.
- Setting ALEMBIC_FOUND to FALSE.\n"
)
if (${ALEMBIC_FOUND})
    if (${IMATH_FOUND})
        target_link_libraries("Alembic" INTERFACE ${THIRDPARTY_IMATH})
    else()
        # Disable Alembic
        set(ALEMBIC_REQUIRED OFF CACHE INTERNAL "")
        set(ALEMBIC_FOUND FALSE CACHE INTERNAL "")
        # Ovewrite default message
        set("ALEMBIC_FOUND_MESSAGE" "ALEMBIC  (ALEMBIC_REQUIRED set to OFF due to missing dependencies)" )
        # print warning
        message(WARNING "${warn_message_libs}")
    endif()
endif()
