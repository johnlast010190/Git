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
    (c) 2019-2024 Engys Ltd.

[----------------------------------------------------------------------------]]

set(TARGET_NAME "dynamicFvMesh")

add_helyx_library(${TARGET_NAME} SHARED
./dynamicFvMesh/dynamicFvMesh.C
./dynamicFvMesh/dynamicFvMeshNew.C
./staticFvMesh/staticFvMesh.C
./dynamicMotionSolverFvMesh/dynamicMotionSolverFvMesh.C
./dynamicMultiMotionSolverFvMesh/dynamicMultiMotionSolverFvMesh.C
./dynamicInkJetFvMesh/dynamicInkJetFvMesh.C
./dynamicFreeSurfaceFvMesh/dynamicFreeSurfaceFvMesh.C
./dynamicRefineFvMesh/dynamicRefineFvMesh.C
./dynamicMultiCriterionRefineFvMesh/dynamicMultiCriterionRefineFvMesh.C
./dynamicSwitchingFvMesh/dynamicSwitchingFvMesh.C
./multiSolidBodyMotionMultiCriterionDynamicRefineFvMesh/multiSolidBodyMotionMultiCriterionDynamicRefineFvMesh.C
./dynamicMotionSolverFvMesh/dynamicMotionSolverFvMesh.C
./dynamicMotionSolverFvMeshAMI/dynamicMotionSolverFvMeshAMI.C
./patchSwitchFvMesh/patchSwitchFvMesh.C
./patchSwitchFvMesh/patchSwitches/patchSwitch/patchSwitch.C
./patchSwitchFvMesh/patchSwitches/patchSwitch/patchSwitchNew.C
./patchSwitchFvMesh/patchSwitches/trapDoor/trapDoor.C
./patchSwitchFvMesh/patchSwitches/fireSensorDoor/fireSensorDoor.C
./patchSwitchFvMesh/patchSwitches/tableSwitch/tableSwitch.C
./patchSwitchFvMesh/patchSwitches/openPatch/openPatch.C
./sliderMotionFvMesh/multiSliderMotionFvMesh.C
./dynamicGIBFvMesh/dynamicGIBFvMesh/dynamicGIBFvMesh.C
./dynamicGIBFvMesh/dynamicGIBFvMesh/dynamicGIBFvMeshBoundary.C
./dynamicGIBFvMesh/dynamicGIBFvMesh/dynamicGIBFvMeshDiagnostics.C
./dynamicGIBFvMesh/dynamicGIBFvMesh/dynamicGIBFvMeshFlipping.C
./dynamicGIBFvMesh/solidBodyGIBFvMesh/solidBodyGIBFvMesh.C
./dynamicGIBFvMesh/deformingBodyGIBFvMesh/deformingBodyGIBFvMesh.C
./dynamicGIBFvMesh/adjointGIBFvMesh/adjointGIBFvMesh.C
./dynamicGIBFvMesh/adjointGIBFvMesh/adjointGIBFvMeshWalkAlgorithm.C
./dynamicGIBFvMesh/movingGIBTools/mapGIB/mapGIB.C
./dynamicGIBFvMesh/movingGIBTools/surfaceAreaSmoothing/surfaceAreaSmoothing.C
./dynamicGIBFvMesh/movingGIBTools/parallelIntersectionData/parallelIntersectionData.C
./dynamicGIBFvMesh/movingGIBTools/constraintBoundaryPoints/constraintBoundaryPoints.C
./dynamicGIBFvMesh/deformingBodyMotionFunction/deformingBodyMotionFunction.C
./dynamicGIBFvMesh/deformingBodyMotionFunction/deformingBodyMotionFunctionNew.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        ${HELYX_PROJECT_DIR}/src/transportModels
        ${HELYX_PROJECT_DIR}/src/transportModels/incompressible/singlePhaseTransportModel
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE meshTools
    PRIVATE dynamicMesh
    PRIVATE finiteVolume
    PRIVATE finiteArea
    PRIVATE OpenFOAM
    PRIVATE incompressibleTransportModels
    PRIVATE turbulenceModels
    PRIVATE surfMesh
    PRIVATE decompositionMethods
    PRIVATE distributed
    PRIVATE ${PTSCOTCHDECOMP}
)
