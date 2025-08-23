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

set(TARGET_NAME "snappyHexMesh")

# Fix issues with FPE for tree vectorisation run in SP on AMD
# Issues in meshing on AMD when using tree vectorize and SP
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND "$ENV{HELYX_PRECISION_OPTION}" MATCHES "SP")
    add_compile_options(-fno-tree-vectorize)
endif()

add_helyx_library(${TARGET_NAME} SHARED
# For the unity build, this file must be included before snappyLayerDriver.C to
# ensure that templates are instantiated in the right order.
./externalDisplacementMeshMover/zeroFixedValue/zeroFixedValuePointPatchFields.C

# snappyRefineDriver.C includes OpenVDB and some math functions (sin, sqrt, etc)
# clash with the OpenFOAM implementation
./helyxMeshGenerator/helyxMeshGenerator.C
./snappyHexMeshDriver/snappyRefineDriver.C
./snappyHexMeshDriver/snappyLayerDriver.C
./snappyHexMeshDriver/snappySnapDriver.C
./snappyHexMeshDriver/snappySnapDriverFeature.C
./snappyHexMeshDriver/snappyVoxelMeshDriver.C
./snappyHexMeshDriver/layerParameters/layerParameters.C
./snappyHexMeshDriver/refinementParameters/refinementParameters.C
./snappyHexMeshDriver/snapParameters/snapParameters.C
./meshRefinement/meshRefinementBaffles.C
./meshRefinement/meshRefinement.C
./meshRefinement/meshRefinementMerge.C
./meshRefinement/meshRefinementProblemCells.C
./meshRefinement/meshRefinementRefine.C
./meshRefinement/meshRefinementGapRefine.C
./meshRefinement/patchFaceOrientation.C
./meshRefinement/fieldData/fieldData.C
./meshRefinement/meshRefinementSplit.C
./meshRefinement/meshRefinementWrapper.C
./snappyHexMeshDriver/pointData/pointData.C
./snappyHexMeshDriver/keepData/keepData.C
./refinementFeatures/manifoldFeatures.C
./refinementFeatures/refinementFeatures.C
./refinementSurfaces/surfaceZonesInfo.C
./refinementSurfaces/refinementSurfaces.C
./shellSurfaces/shellSurfaces.C
./trackedParticle/trackedParticle.C
./trackedParticle/trackedParticleCloud.C
./featureLineSnapping/featureLinePrep.C
./featureLineSnapping/featureLineSnapper.C
./hexReport/hexReport.C
./addHexInternalLayers/addHexInternalLayers.C
./addHexMeshLayers/addHexMeshLayer.C
./addHexMeshLayers/addMultiLayers.C
./autoHexDualiser/autoHexDualiser.C
./autoHexDualiser/extendedEdgeCells.C
./autoHexDualiser/extendedEdgeCellsIO.C
./autoHexDualiser/extendedPointEdgeFaces.C
./autoHexDualiser/extendedPointEdgeFacesIO.C
./autoBlockMesh/autoBlockMesh.C
./autoFanExtrude/autoFanExtrude.C
./autoExtrude/autoExtrude.C
./autoLayerCellsMerge/autoLayerCellsMerge.C
./layerManipulate/layerManipulate.C
./externalDisplacementMeshMover/displacementMeshMoverMotionSolver.C
./externalDisplacementMeshMover/externalDisplacementMeshMover.C
./externalDisplacementMeshMover/medialAxisMeshMover.C
./externalDisplacementMeshMover/displacementMotionSolverMeshMover.C
./externalDisplacementMeshMover/fieldSmoother/fieldSmoother.C
./anisoRefiner/anisoRefiner.C
./autoSplitCells/autoSplitCells.C
./autoOptimize/autoOptimize.C
./autoOptimize/autoOptimizeNew.C
./autoOptimize/noOptimize/noOptimize.C
./autoOptimize/cfMeshOptimize/cfMeshOptimize.C
./autoOptimize/helyxOptimize/helyxOptimize.C
./edgeClassification/edgeClassification.C
./errorCellRemoval/errorCellRemoval.C
./meshControl/meshControl.C
./hessianMeshOptimization/leastSquaresCurveFit/leastSquaresCurveFit.C
./hessianMeshOptimization/meshQualityStatus/meshQualityStatus.C
./hessianMeshOptimization/meshOptimizationMethod/meshOptimizationMethod.C
./hessianMeshOptimization/meshOptimizationMethod/meshOptimizationMethodNew.C
./hessianMeshOptimization/meshOptimizationMethod/activeSet/activeSet.C
./hessianMeshOptimization/meshOptimizationMethod/fullMeshOptimization/fullMeshOptimization.C
./hessianMeshOptimization/meshOptimizationMethod/qualityBasedOptimization/qualityBasedOptimization.C
./hessianMeshOptimization/meshGeometry/meshGeometry.C
./hessianMeshOptimization/meshMetric/meshMetric.C
./hessianMeshOptimization/layerCellHandling/stencilBuilder/stencilBuilder.C
./hessianMeshOptimization/layerCellHandling/stencilBuilder/stencilBuilderNew.C
./hessianMeshOptimization/layerCellHandling/stencilBuilder/isotropy/isotropy.C
./hessianMeshOptimization/layerCellHandling/stencilBuilder/layerAnisotropy/layerAnisotropy.C
./hessianMeshOptimization/layerCellHandling/layerCellDetection.C
./hessianMeshOptimization/layerCellHandling/cellMetric.C
./hessianMeshOptimization/layerCellHandling/layerPointStack/layerPointStack.C
./hessianMeshOptimization/layerCellHandling/layerPointStack/layerPointStackIO.C
./hessianMeshOptimization/smoother/generalMeshSmoother.C
./hessianMeshOptimization/meshOptimization/hessianMeshOptimization.C
./orientFaceZone/orientFaceZone.C
./boundaryLayerRefinement/boundaryLayerRefinement.C
)


helyx_link_libraries(${TARGET_NAME}
    PRIVATE meshLibrary
    PRIVATE dynamicMesh
    PRIVATE finiteVolume
    PRIVATE lagrangian
    PRIVATE meshTools
    PRIVATE surfMesh
    PRIVATE fvMotionSolvers
    PRIVATE distributed
    PRIVATE sampling
    PRIVATE overset
    PRIVATE fileFormats
    PRIVATE decompose
    PRIVATE decompositionMethods
    PUBLIC ${SCOTCHDECOMP}
    PUBLIC ${PTSCOTCHDECOMP}
    PUBLIC ${PARHIPDECOMP}
    PRIVATE renumberMethods
    PRIVATE blockMesh
)

if("${OPENVDB_FOUND}")
    add_compile_definitions(HELYX_USE_OPENVDB)

    helyx_link_libraries(${TARGET_NAME}
        #PRIVATE ${THIRDPARTY_TBB}     # It is already linked to OpenVDB
        #PRIVATE ${THIRDPARTY_BOOST}   # It is already linked to OpenVDB
        PRIVATE ${THIRDPARTY_OPENVDB}  # Do we need this? helyxVDB already links against OpenVDB
        PRIVATE helyxVDB
    )
endif()
