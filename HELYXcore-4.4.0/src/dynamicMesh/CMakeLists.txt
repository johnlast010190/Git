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
    (c) 2019-2023 Engys Ltd.

[----------------------------------------------------------------------------]]

set(TARGET_NAME "dynamicMesh")

add_helyx_library(${TARGET_NAME} SHARED
./attachDetach/attachDetach.C
./attachDetach/attachInterface.C
./attachDetach/detachInterface.C
./attachDetach/attachDetachPointMatchMap.C
./layerAdditionRemoval/layerAdditionRemoval.C
./layerAdditionRemoval/setLayerPairing.C
./layerAdditionRemoval/addCellLayer.C
./layerAdditionRemoval/removeCellLayer.C
./slidingInterface/enrichedPatch/enrichedPatch.C
./slidingInterface/enrichedPatch/enrichedPatchPointMap.C
./slidingInterface/enrichedPatch/enrichedPatchFaces.C
./slidingInterface/enrichedPatch/enrichedPatchPointPoints.C
./slidingInterface/enrichedPatch/enrichedPatchCutFaces.C
./slidingInterface/enrichedPatch/enrichedPatchMasterPoints.C
./slidingInterface/slidingInterface.C
./slidingInterface/slidingInterfaceProjectPoints.C
./slidingInterface/coupleSlidingInterface.C
./slidingInterface/slidingInterfaceAttachedAddressing.C
./slidingInterface/slidingInterfaceClearCouple.C
./slidingInterface/decoupleSlidingInterface.C
./perfectInterface/perfectInterface.C
./boundaryMesh/boundaryMesh.C
./boundaryPatch/boundaryPatch.C
./setUpdater/setUpdater.C
./meshCut/meshModifiers/boundaryCutter/boundaryCutter.C
./meshCut/meshModifiers/meshCutter/meshCutter.C
./meshCut/meshModifiers/meshCutAndRemove/meshCutAndRemove.C
./meshCut/meshModifiers/undoableMeshCutter/undoableMeshCutter.C
./meshCut/meshModifiers/refinementIterator/refinementIterator.C
./meshCut/meshModifiers/multiDirRefinement/multiDirRefinement.C
./meshCut/cellLooper/cellLooper.C
./meshCut/cellLooper/topoCellLooper.C
./meshCut/cellLooper/geomCellLooper.C
./meshCut/cellLooper/hexCellLooper.C
./meshCut/directions/directions.C
./meshCut/directions/directionInfo/directionInfo.C
./meshCut/edgeVertex/edgeVertex.C
./meshCut/cellCuts/cellCuts.C
./meshCut/splitCell/splitCell.C
./meshCut/refineCell/refineCell.C
./meshCut/wallLayerCells/wallLayerCells.C
./meshCut/wallLayerCells/wallNormalInfo/wallNormalInfo.C
./fvMeshAdder/fvMeshAdder.C
./fvMeshDistribute/fvMeshDistribute.C
./meshDistribute/meshAdder.C
./meshDistribute/meshDistribute.C
./fvMeshDistribute/IOmapDistributePolyMesh.C
./polyMeshAdder/faceCoupleInfo.C
./polyMeshAdder/polyMeshAdder.C
./fvMeshTools/fvMeshTools.C
./fvMeshSubset/fvMeshSubset.C
./fvMeshSubsetProxy/fvMeshSubsetProxy.C
./meshSubsetHelper/meshSubsetHelper.C
./motionSmoother/motionSmoother.C
./motionSmoother/motionSmootherAlgo.C
./motionSmoother/motionSmootherAlgoCheck.C
./motionSmoother/motionSmootherData.C
./motionSmoother/polyMeshGeometry/polyMeshGeometry.C
./motionSmoother/polyMeshGeometry/meshStatistics.C
./motionSmoother/badQualityToCell/badQualityToCell.C
./motionSmoother/badQualityToFace/badQualityToFace.C
./motionSolvers/motionSolver/motionSolver.C
./motionSolvers/displacement/points/pointsMotionSolver.C
./motionSolvers/displacement/displacement/displacementMotionSolver.C
./motionSolvers/displacement/interpolation/displacementInterpolationMotionSolver.C
./motionSolvers/displacement/layeredSolver/displacementLayeredMotionMotionSolver.C
./motionSolvers/displacement/layeredSolver/pointEdgeStructuredWalk.C
./motionSolvers/componentDisplacement/componentDisplacementMotionSolver.C
./motionSolvers/velocity/velocityMotionSolver.C
./motionSolvers/velocity/velocityDisplacement/velocityDisplacementMotionSolver.C
./motionSolvers/componentVelocity/componentVelocityMotionSolver.C
./motionSolvers/displacement/solidBody/solidBodyMotionSolver.C
./motionSolvers/displacement/solidBody/multiSolidBodyMotionSolver.C
./motionSolvers/displacement/solidBody/multiSolidBodyGIBMotionSolver.C
./motionSolvers/displacement/solidBody/pointPatchFields/derived/solidBodyMotionDisplacement/solidBodyMotionDisplacementPointPatchVectorField.C
./motionSolvers/modes/modesMotionSolver.C
./motionSolvers/sinus/sinusMotionSolver.C
./motionSolvers/motionSolverList/motionSolverList.C
./createShellMesh/createShellMesh.C
./extrudePatchMesh/extrudePatchMesh.C
./polyMeshFilter/polyMeshFilterSettings.C
./polyMeshFilter/polyMeshFilter.C
./pointPatchDist/externalPointEdgePoint.C
./pointPatchDist/pointPatchDist.C
./polyTopoChange/edgeCollapser/edgeCollapser.C
./polyTopoChange/faceCollapser/faceCollapser.C
./polyTopoChange/pointEdgeCollapse/pointEdgeCollapse.C
./polyTopoChange/removeCells/removeCells.C
./polyTopoChange/removeFaces/removeFaces.C
./polyTopoChange/removePoints/removePoints.C
./polyTopoChange/refinementData/refinementData.C
./polyTopoChange/refinementDistanceData/refinementDistanceData.C
./polyTopoChange/combineFaces/combineFaces.C
./polyTopoChange/duplicatePoints/duplicatePoints.C
./polyTopoChange/tetDecomposer/tetDecomposer.C
./polyTopoChange/hexRef8/hexRef8.C
./polyTopoChange/hexRef8/hexRef8Data.C
./polyTopoChange/hexRef8/hexRef8DataList.C
./polyTopoChange/hexRef8/refinementHistory.C
)



helyx_link_libraries(${TARGET_NAME}
    PRIVATE finiteVolume
    PRIVATE surfMesh
    PRIVATE extrudeModel
    PUBLIC meshTools
)

if("${OPENVDB_FOUND}")
    helyx_link_libraries(${TARGET_NAME}
        #PRIVATE ${THIRDPARTY_TBB}   # OpenVDB is already linking against
        #PRIVATE ${THIRDPARTY_BOOST} # OpenVDB is already linking against
        PRIVATE ${THIRDPARTY_OPENVDB}
    )
endif()
