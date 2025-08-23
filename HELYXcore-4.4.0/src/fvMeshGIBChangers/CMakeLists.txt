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
    (c) 2024 Engys Ltd.

[----------------------------------------------------------------------------]]

set(TARGET_NAME "fvMeshGIBChangers")

add_helyx_library(${TARGET_NAME} SHARED
./base/fvMeshGIBChangersBase.C
./base/fvMeshGIBChangersBaseBoundary.C
./base/fvMeshGIBChangersBaseDiagnostics.C
./base/fvMeshGIBChangersBaseFlipping.C
./GIBTools/GIBAreaSmoothing/GIBAreaSmoothing.C
./GIBTools/GIBIntersectionData/GIBIntersectionData.C
./GIBTools/GIBMapping/GIBMapping.C
./deformingBody/fvMeshGIBChangersDeformingBody.C
./deformingBody/motionFunctions/motionFunction/deformingBodyMotionFunction.C
./deformingBody/motionFunctions/motionFunction/deformingBodyMotionFunctionNew.C
./deformingBody/motionFunctions/adjointMotion/adjointMotion.C
./deformingBody/motionFunctions/burningMotion/burningMotion.C
./deformingBody/motionFunctions/staticMotion/staticMotion.C
./adjoint/fvMeshGIBChangersAdjoint.C
./adjoint/fvMeshGIBChangersAdjointWalkAlgorithm.C
./patchSwitch/fvMeshGIBChangersPatchSwitch.C
./patchSwitch/switches/GIBSwitch/GIBSwitch.C
./patchSwitch/switches/GIBSwitch/GIBSwitchNew.C
./patchSwitch/switches/fireSensorDoorGIBSwitch/fireSensorDoorGIBSwitch.C
./patchSwitch/switches/openPatch/openPatch.C
./patchSwitch/switches/tableSwitch/tableSwitch.C
./patchSwitch/switches/trapDoor/trapDoor.C
./solidBody/fvMeshGIBChangersSolidBody.C
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE finiteVolume
    PRIVATE dynamicMesh
    PRIVATE finiteArea
    PRIVATE fvMeshStitchers
)

if (${swak4Foam_FOUND})
    helyx_link_libraries(${TARGET_NAME} PRIVATE swakSourceFields)
    helyx_link_libraries_with_generated_headers(${TARGET_NAME} PRIVATE swak4FoamParsers)
else()
    message(STATUS "Disabling expression support in burningMotion (swak4Foam not found)")
    target_compile_definitions(${TARGET_NAME}_obj PRIVATE HELYX_SWAK4FOAM_DISABLED)
endif ()
