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
Custom findModule for OpenCASCADE

[----------------------------------------------------------------------------]]


# N.B.  OpenCASCADE headers are not in standard include dir (*/include/),
# instead, OpenCASCADE project assumes that headers base dir is */include/opencascade/.
set(opencascade_libs
"TKernel;\
TKLCAF;\
TKCDF;\
TKService;\
TKXmlL;\
TKBinL;\
TKXml;\
TKVCAF;\
TKXmlXCAF;\
TKRWMesh;\
TKMath;\
TKBinXCAF;\
TKBin;\
TKV3d;\
TKXCAF;\
TKHLR;\
TKCAF;\
TKBRep;\
TKG3d;\
TKTopAlgo;\
TKXSBase;\
TKMesh;\
TKShHealing;\
TKGeomAlgo;\
TKGeomBase;\
TKG2d;\
TKBool;\
TKBO;\
TKDEIGES;\
TKDESTEP;\
TKDESTL;\
TKPrim")

helyx_find_thirdparty_package(OPENCASCADE
    "${opencascade_libs}"
    "STEPControl_Reader.hxx")

# If OpenCASCADE compilation is disabled, simply return from this script.
# There is no need to evaluate dependencies
if(OPENCASCADE_REQUIRED STREQUAL "OFF" OR NOT OPENCASCADE_FOUND)
    return()
endif()

## N.B.  Could add extra checks and detection here to not have TBB
## on OPTIONAL_THIRDPARTY_LIBRARIES:
## Check combinations of OPENCASCADE_REQUIRED and dependency_REQUIRED
## then set ${OPENCASCADE|dependency}_FOUND or ${OPENCASCADE|dependency}_REQUIRED accordingly
#helyx_check_module_mandatory_dependency(OPENCASCADE TBB)

# OpenCascade optionally depends on TBB.
# Link TBB if found, otherwise, send an warning message.
# At this point, all dependencies should have been evaluated.
STRING(CONCAT warn_message_libs
"OPENCASCADE libraries found, but optional dependencies are missing:
    TBB_FOUND: ${TBB_FOUND}
If you are willing to use TBB, review the dependencies setup in the settings file and reconfigure.
To enable 'TBB', check if:
    TBB is built in ThirdParty or available in the system,
    TBB_ARCH_PATH is properly set in the userSettings file,
    TBB_REQUIRED is set to ON|AUTO.\n"
)
if (${OPENCASCADE_FOUND})
    if(${TBB_FOUND})
        foreach(lib IN LISTS opencascade_libs)
            target_link_libraries("${lib}" INTERFACE ${THIRDPARTY_TBB})
        endforeach()
    else()
        ## Disable OpenCascade
        #set(OPENCASCADE_REQUIRED OFF CACHE INTERNAL "")
        #set(OPENCASCADE_FOUND FALSE CACHE INTERNAL "")
        ## Ovewrite default message
        #set("OPENCASCADE_FOUND_MESSAGE" "OPENCASCADE  (OPENCASCADE_REQUIRED set to OFF due to missing dependencies)" )
        message(WARNING "${warn_message_libs}")
    endif()
endif()
