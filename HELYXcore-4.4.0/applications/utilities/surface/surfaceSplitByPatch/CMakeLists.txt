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

[----------------------------------------------------------------------------]]

set(TARGET_NAME "surfaceSplitByPatch-exe")

add_helyx_executable(${TARGET_NAME}
./surfaceSplitByPatch.C
)



helyx_link_libraries(${TARGET_NAME}
    PRIVATE surfMesh
    PRIVATE OpenFOAM
    PRIVATE decompose
)

# Runtime library loading for Windows
if(${HELYX_SYSTEM_NAME} STREQUAL "MSwindows")
    #get_target_property(LINK_LIBRARIES_STRING ${TARGET_NAME} LINK_LIBRARIES)
    #string(CONCAT LIBS_TO_LOAD "\"${LINK_LIBRARIES_STRING}\"")
    #string (REPLACE ";" "," LIBS_TO_LOAD "${LIBS_TO_LOAD}")
    #string (REGEX REPLACE "\/.*\/" "" LIBS_TO_LOAD "${LIBS_TO_LOAD}")
    #target_compile_definitions(${TARGET_NAME} PRIVATE "LIBS_TO_LOAD=${LIBS_TO_LOAD}")

    set_target_properties(${TARGET_NAME} PROPERTIES OUTPUT_NAME surfaceSplitByPtch)
    file(WRITE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/surfaceSplitByPatch.bat "@surfaceSplitByPtch.exe %*")
endif()
