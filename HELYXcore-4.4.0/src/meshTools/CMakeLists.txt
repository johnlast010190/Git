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
    (c) 2019-2025 Engys Ltd.

[----------------------------------------------------------------------------]]

set(TARGET_NAME "meshTools")

add_helyx_library(${TARGET_NAME} SHARED
./edgeMesh/edgeMesh.C
./edgeMesh/edgeMeshIO.C
./edgeMesh/edgeMeshNew.C
./edgeMesh/edgeMeshFormats/edgeMeshFormatsCore.C
./edgeMesh/edgeMeshFormats/edgeMesh/edgeMeshFormat.C
./edgeMesh/edgeMeshFormats/edgeMesh/edgeMeshFormatRunTime.C
./edgeMesh/edgeMeshFormats/extendedFeatureEdgeMesh/extendedFeatureEdgeMeshFormat.C
./edgeMesh/edgeMeshFormats/extendedFeatureEdgeMesh/extendedFeatureEdgeMeshFormatRunTime.C
./edgeMesh/edgeMeshFormats/nas/NASedgeFormat.C
./edgeMesh/edgeMeshFormats/nas/NASedgeFormatRunTime.C
./edgeMesh/edgeMeshFormats/obj/OBJedgeFormat.C
./edgeMesh/edgeMeshFormats/obj/OBJedgeFormatRunTime.C
./edgeMesh/edgeMeshFormats/starcd/STARCDedgeFormat.C
./edgeMesh/edgeMeshFormats/starcd/STARCDedgeFormatRunTime.C
./edgeMesh/edgeMeshFormats/vtk/VTKedgeFormat.C
./edgeMesh/edgeMeshFormats/vtk/VTKedgeFormatRunTime.C
./edgeMesh/edgeMeshTools/edgeMeshTools.C
./edgeMesh/edgeMeshTools/edgeMeshFeatureProximity.C
./edgeMesh/featureEdgeMesh/featureEdgeMesh.C
./edgeMesh/extendedEdgeMesh/extendedEdgeMesh.C
./edgeMesh/extendedEdgeMesh/extendedEdgeMeshNew.C
./edgeMesh/extendedEdgeMesh/extendedEdgeMeshFormats/extendedEdgeMeshFormat/extendedEdgeMeshFormat.C
./edgeMesh/extendedEdgeMesh/extendedEdgeMeshFormats/extendedEdgeMeshFormat/extendedEdgeMeshFormatRunTime.C
./edgeMesh/extendedEdgeMesh/extendedFeatureEdgeMesh/extendedFeatureEdgeMesh.C
./cellClassification/cellClassification.C
./cellClassification/cellInfo.C
./cellQuality/cellQuality.C
./patchDist/patchDistFuncs/patchDistFuncs.C
./patchDist/wallPoint/wallPoint.C
./cellFeatures/cellFeatures.C
./polyMeshZipUpCells/polyMeshZipUpCells.C
./primitiveMeshGeometry/primitiveMeshGeometry.C
./meshSearch/meshSearch.C
./meshSearch/meshSearchFACE_CENTRE_TRISMeshObject.C
./meshSearch/meshSearchMeshObject.C
./meshTools/meshTools.C
./algorithms/PointEdgeWave/PointEdgeWaveName.C
./algorithms/PointEdgeWave/pointEdgePoint.C
./algorithms/PatchEdgeFaceWave/PatchEdgeFaceWaveName.C
./algorithms/PatchEdgeFaceWave/patchEdgeFaceInfo.C
./algorithms/PatchEdgeFaceWave/patchEdgeFaceBoolInfo.C
./algorithms/PatchEdgeFaceWave/patchPatchDist.C
./algorithms/PatchEdgeFaceWave/patchEdgeFaceRegion.C
./algorithms/PatchEdgeFaceWave/patchEdgeFaceRegions.C
./algorithms/FaceCellWave/FaceCellWaveName.C
./regionSplit/regionSplit.C
./regionSplit/localPointRegion.C
./regionSplit/regionSplitPatch.C
./regionSplit2D/regionSplit2D.C
./indexedOctree/treeDataEdge.C
./indexedOctree/treeDataFace.C
./indexedOctree/treeDataPoint.C
./indexedOctree/treeDataPrimitivePatchName.C
./indexedOctree/treeDataTriSurface.C
./searchableSurfaces/closedTriSurfaceMesh/closedTriSurfaceMesh.C
./searchableSurfaces/searchableBox/searchableBox.C
./searchableSurfaces/searchableCone/searchableCone.C
./searchableSurfaces/searchableCylinder/searchableCylinder.C
./searchableSurfaces/searchableRing/searchableRing.C
./searchableSurfaces/searchableDisk/searchableDisk.C
./searchableSurfaces/searchableExtrudedCircle/searchableExtrudedCircle.C
./searchableSurfaces/searchablePlane/searchablePlane.C
./searchableSurfaces/searchablePlate/searchablePlate.C
./searchableSurfaces/searchableRotatedBox/searchableRotatedBox.C
./searchableSurfaces/searchableSphere/searchableSphere.C
./searchableSurfaces/searchableSurface/searchableSurface.C
./searchableSurfaces/searchableSurfaceCollection/searchableSurfaceCollection.C
./searchableSurfaces/searchableSurfaces/searchableSurfaces.C
./searchableSurfaces/searchableSurfacesQueries/searchableSurfacesQueries.C
./searchableSurfaces/searchableSurfaceWithGaps/searchableSurfaceWithGaps.C
./searchableSurfaces/subTriSurfaceMesh/subTriSurfaceMesh.C
./searchableSurfaces/triSurfaceMesh/triSurfaceMesh.C
./sets/topoSets/cellSet.C
./sets/topoSets/topoSet.C
./sets/topoSets/faceSet.C
./sets/topoSets/pointSet.C
./sets/topoSets/cellZoneSet.C
./sets/topoSets/faceZoneSet.C
./sets/topoSets/pointZoneSet.C
./sets/topoSetSource/topoSetSource.C
./sets/cellSources/faceToCell/faceToCell.C
./sets/cellSources/fieldToCell/fieldToCell.C
./sets/cellSources/pointToCell/pointToCell.C
./sets/cellSources/shapeToCell/shapeToCell.C
./sets/cellSources/boxToCell/boxToCell.C
./sets/cellSources/regionToCell/regionToCell.C
./sets/cellSources/rotatedBoxToCell/rotatedBoxToCell.C
./sets/cellSources/labelToCell/labelToCell.C
./sets/cellSources/surfaceToCell/surfaceToCell.C
./sets/cellSources/cellToCell/cellToCell.C
./sets/cellSources/nearestToCell/nearestToCell.C
./sets/cellSources/nbrToCell/nbrToCell.C
./sets/cellSources/zoneToCell/zoneToCell.C
./sets/cellSources/sphereToCell/sphereToCell.C
./sets/cellSources/cylinderToCell/cylinderToCell.C
./sets/cellSources/faceZoneToCell/faceZoneToCell.C
./sets/cellSources/cylinderAnnulusToCell/cylinderAnnulusToCell.C
./sets/cellSources/patchToCell/patchToCell.C
./sets/cellSources/targetVolumeToCell/targetVolumeToCell.C
./sets/faceSources/faceToFace/faceToFace.C
./sets/faceSources/labelToFace/labelToFace.C
./sets/faceSources/cellToFace/cellToFace.C
./sets/faceSources/normalToFace/normalToFace.C
./sets/faceSources/pointToFace/pointToFace.C
./sets/faceSources/patchToFace/patchToFace.C
./sets/faceSources/boundaryToFace/boundaryToFace.C
./sets/faceSources/zoneToFace/zoneToFace.C
./sets/faceSources/boxToFace/boxToFace.C
./sets/faceSources/regionToFace/regionToFace.C
./sets/faceSources/cylinderToFace/cylinderToFace.C
./sets/pointSources/labelToPoint/labelToPoint.C
./sets/pointSources/pointToPoint/pointToPoint.C
./sets/pointSources/cellToPoint/cellToPoint.C
./sets/pointSources/faceToPoint/faceToPoint.C
./sets/pointSources/boxToPoint/boxToPoint.C
./sets/pointSources/surfaceToPoint/surfaceToPoint.C
./sets/pointSources/zoneToPoint/zoneToPoint.C
./sets/pointSources/nearestToPoint/nearestToPoint.C
./sets/faceZoneSources/faceZoneToFaceZone/faceZoneToFaceZone.C
./sets/faceZoneSources/setsToFaceZone/setsToFaceZone.C
./sets/faceZoneSources/setToFaceZone/setToFaceZone.C
./sets/faceZoneSources/setAndNormalToFaceZone/setAndNormalToFaceZone.C
./sets/faceZoneSources/searchableSurfaceToFaceZone/searchableSurfaceToFaceZone.C
./sets/faceZoneSources/planeToFaceZone/planeToFaceZone.C
./sets/cellZoneSources/setToCellZone/setToCellZone.C
./sets/pointZoneSources/setToPointZone/setToPointZone.C
./momentOfInertia/momentOfInertia.C
./centreOfMass/centreOfMass.C
./surfaceSets/surfaceSets.C
./triSurface/orientedSurface/orientedSurface.C
./triSurface/surfaceLocation/surfaceLocation.C
./triSurface/booleanOps/surfaceIntersection/surfaceIntersection.C
./triSurface/booleanOps/surfaceIntersection/surfaceIntersectionFuncs.C
./triSurface/booleanOps/surfaceIntersection/edgeIntersections.C
./triSurface/booleanOps/booleanSurface/booleanSurface.C
./triSurface/booleanOps/intersectedSurface/intersectedSurface.C
./triSurface/booleanOps/intersectedSurface/edgeSurface.C
./triSurface/faceTriangulation/faceTriangulation.C
./triSurface/triSurfaceSearch/triSurfaceSearch.C
./triSurface/triSurfaceSearch/triSurfaceRegionSearch.C
./triSurface/triangleFuncs/triangleFuncs.C
./triSurface/surfaceFeatures/surfaceFeatures.C
./triSurface/triSurfaceLoader/triSurfaceLoader.C
./triSurface/triSurfaceTools/triSurfaceTools.C
./triSurface/triSurfaceTools/triSurfaceCloseness.C
./triSurface/triSurfaceTools/triSurfaceCurvature.C
./triSurface/triSurfaceTools/geompack/geompack.C
./triSurface/triSurfaceTools/pointToPointPlanarInterpolation.C
./twoDPointCorrector/twoDPointCorrector.C
./meshQualityMetrics/meshQualityMetrics.C
./meshQualityMetrics/meshQualityMetricsNew.C
./AMIInterpolation//AMIInterpolation/AMIInterpolation.C
./AMIInterpolation//AMIInterpolation/AMIInterpolationParallelOps.C
./AMIInterpolation/faceAreaIntersect/faceAreaIntersect.C
./AMIInterpolation/GAMG/interfaces/cyclicAMIGAMGInterface/cyclicAMIGAMGInterface.C
./AMIInterpolation/GAMG/interfaceFields/cyclicAMIGAMGInterfaceField/cyclicAMIGAMGInterfaceField.C
./AMIInterpolation/AMIInterpolation/AMIMethod/AMIMethod/AMIMethod.C
./AMIInterpolation/AMIInterpolation/AMIMethod/AMIMethod/AMIMethodNew.C
./AMIInterpolation/AMIInterpolation/AMIMethod/directAMI/directAMI.C
./AMIInterpolation/AMIInterpolation/AMIMethod/faceAreaWeightAMI/faceAreaWeightAMI.C
./AMIInterpolation/AMIInterpolation/AMIMethod/mapNearestAMI/mapNearestAMI.C
./AMIInterpolation/AMIInterpolation/AMIMethod/partialFaceAreaWeightAMI/partialFaceAreaWeightAMI.C
./AMIInterpolation/patches/cyclicAMI/cyclicAMILduInterfaceField/cyclicAMILduInterface.C
./AMIInterpolation/patches/cyclicAMI/cyclicAMILduInterfaceField/cyclicAMILduInterfaceField.C
./AMIInterpolation/patches/cyclicAMI/cyclicAMIPolyPatch/cyclicAMIPolyPatch.C
./AMIInterpolation/patches/cyclicAMI/cyclicAMIPolyPatch/cyclicAMIPolyPatchTopologyChange.C
./AMIInterpolation/patches/cyclicAMI/cyclicAMIPointPatch/cyclicAMIPointPatch.C
./AMIInterpolation/patches/cyclicAMI/cyclicAMIPointPatchField/cyclicAMIPointPatchFields.C
./AMIInterpolation/patches/cyclicPeriodicAMI/cyclicPeriodicAMIPolyPatch/cyclicPeriodicAMIPolyPatch.C
./AMIInterpolation/patches/discreteMixingPlane/discreteMixingPlanePolyPatch/discreteMixingPlanePolyPatch.C
./mappedPatches/mappedPolyPatch/mappedPatchBase.C
./mappedPatches/mappedPolyPatch/mappedPolyPatch.C
./mappedPatches/mappedPolyPatch/mappedWallPolyPatch.C
./mappedPatches/mappedPointPatch/mappedPointPatch.C
./mappedPatches/mappedPointPatch/mappedWallPointPatch.C
./mappedPatches/mappedLduInterface/mappedLduInterface.C
./mappedPatches/GAMG/mappedGAMGInterface.C
./mappedPatches/GAMG/mappedGAMGInterfaceField.C
./mappedPatches/GAMG/mappedWallGAMGInterface.C
./mappedPatches/GAMG/mappedWallGAMGInterfaceField.C
./meshStructure/meshStructure.C
./meshStructure/topoDistanceData.C
./meshStructure/pointTopoDistanceData.C
./triIntersect/triIntersect.C
./triIntersect/triIntersectLocationIO.C
./patchToPatch/patchToPatch/patchToPatch.C
./patchToPatch/patchToPatch/patchToPatchParallelOps.C
./patchToPatch/intersection/intersectionPatchToPatch.C
./patchToPatch/inverseDistance/inverseDistancePatchToPatch.C
./patchToPatch/nearby/nearbyPatchToPatch.C
./patchToPatch/nearest/nearestPatchToPatch.C
./patchToPatch/rays/raysPatchToPatch.C
./patchToPatchTools/patchToPatchTools.C
./nonConformal/boundary/nonConformalBoundary.C
./nonConformal/polyPatches/nonConformal/nonConformalPolyPatch.C
./nonConformal/polyPatches/nonConformalCoupled/nonConformalCoupledPolyPatch.C
./nonConformal/polyPatches/nonConformalCyclic/nonConformalCyclicPolyPatch.C
./nonConformal/polyPatches/nonConformalDiscreteMixing/nonConformalDiscreteMixingPolyPatch.C
./nonConformal/polyPatches/nonConformalDiscreteMixing/intersection/discreteMixingIntersection.C
./nonConformal/polyPatches/nonConformalError/nonConformalErrorPolyPatch.C
./nonConformal/polyPatches/nonConformalProcessorCyclic/nonConformalProcessorCyclicPolyPatch.C
./nonConformal/pointPatches/nonConformalCyclic/nonConformalCyclicPointPatch.C
./nonConformal/pointPatches/nonConformalDiscreteMixing/nonConformalDiscreteMixingPointPatch.C
./nonConformal/pointPatches/nonConformalError/nonConformalErrorPointPatch.C
./nonConformal/pointPatches/nonConformalProcessorCyclic/nonConformalProcessorCyclicPointPatch.C
./nonConformal/pointPatchFields/nonConformalCyclic/nonConformalCyclicPointPatchFields.C
./nonConformal/pointPatchFields/nonConformalDiscreteMixing/nonConformalDiscreteMixingPointPatchFields.C
./nonConformal/pointPatchFields/nonConformalError/nonConformalErrorPointPatchFields.C
./nonConformal/pointPatchFields/nonConformalProcessorCyclic/nonConformalProcessorCyclicPointPatchFields.C
./nonConformal/GAMG/nonConformalDiscreteMixing/nonConformalDiscreteMixingGAMGInterface.C
./nonConformal/GAMG/nonConformalDiscreteMixing/nonConformalDiscreteMixingGAMGInterfaceField.C
./output/foamVtkIndPatchWriter.C
./output/foamVtkWriteTopoSet.C
./output/foamVtkWriteFaceSet.C
./output/foamVtkWritePointSet.C
./output/foamVtkWriteCellSetFaces.C
./regionCoupled/patches/regionCoupledPolyPatch/regionCoupledBase.C
./regionCoupled/patches/regionCoupledPolyPatch/regionCoupledPolyPatch.C
./regionCoupled/patches/regionCoupledPolyPatch/regionCoupledWallPolyPatch.C
./regionCoupled/patches/regionCoupledPointPatch/regionCoupledPointPatch.C
./regionCoupled/patches/regionCoupledPointPatch/regionCoupledWallPointPatch.C
./tetOverlapVolume/tetOverlapVolume.C
./vNurbsDeformation/volumetricNurbs/volumetricNurbs.C
./vNurbsDeformation/controlBox/controlBox.C
./vNurbsDeformation/splineBasis/splineBasis.C
./polyTopoChange/polyMeshModifier/polyMeshModifier.C
./polyTopoChange/polyMeshModifier/polyMeshModifierNew.C
./polyTopoChange/polyTopoChange/topoAction/topoActions.C
./polyTopoChange/polyTopoChanger/polyTopoChanger.C
./polyTopoChange/polyTopoChange/polyTopoChange.C
./polyTopoChange/polyTopoChange/addPatchCellLayer.C
./polyTopoChange/attachPolyTopoChanger/attachPolyTopoChanger.C
./polyTopoChange/repatchPolyTopoChanger/repatchPolyTopoChanger.C
./GIB/geometryTransformation/geometryTransformation.C
./GIB/geometryTransformation/geometryTransformationNew.C
./GIB/intersectionData/intersectionData.C
./GIB/twoDPointGIBCorrector/twoDPointGIBCorrector.C
./meshDualiser/meshDualiser.C
)

if("${OPENVDB_FOUND}")
    helyx_link_libraries(${TARGET_NAME}
        PRIVATE ${THIRDPARTY_OPENVDB}
    )
endif()

if("${OPENCASCADE_FOUND}")
    target_sources(${TARGET_NAME}
        PRIVATE ./searchableSurfaces/CADSurface/CADSurface.C
    )
    helyx_link_libraries(${TARGET_NAME}
        PRIVATE ${THIRDPARTY_OPENCASCADE}
    )
endif()

helyx_link_libraries(${TARGET_NAME}
    PUBLIC fileFormats
    PUBLIC surfMesh
    PUBLIC delabella
    PUBLIC OpenFOAM
)

set_source_files_properties(
    # These include file-local functions
    ${CMAKE_CURRENT_LIST_DIR}/coordinate/systems/sphericalCS.C
    ${CMAKE_CURRENT_LIST_DIR}/triSurface/booleanOps/intersectedSurface/edgeSurface.C
    PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON
)
