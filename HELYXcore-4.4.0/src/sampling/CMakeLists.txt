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
    (c) 2019-2022 Engys Ltd.

[----------------------------------------------------------------------------]]

set(TARGET_NAME "sampling")

add_helyx_library(${TARGET_NAME} SHARED
./probes/probes.C
./probes/patchProbes.C
./probes/multiFacePatchProbes.C
./probes/probesGrouping.C
./sampledSet/circle/circleSet.C
./sampledSet/cloud/cloudSet.C
./sampledSet/patchCloud/patchCloudSet.C
./sampledSet/polyLine/polyLineSet.C
./sampledSet/face/faceOnlySet.C
./sampledSet/midPoint/midPointSet.C
./sampledSet/midPointAndFace/midPointAndFaceSet.C
./sampledSet/patchSeed/patchSeedSet.C
./sampledSet/sampledSet/sampledSet.C
./sampledSet/sampledSets/sampledSets.C
./sampledSet/sampledSets/sampledSetsGrouping.C
./sampledSet/triSurfaceMeshPointSet/triSurfaceMeshPointSet.C
./sampledSet/uniform/uniformSet.C
./sampledSet/array/arraySet.C
./sampledSet/shortestPath/shortestPathSet.C
./surface/cuttingPlane/cuttingPlane.C
./surface/cuttingTriSurface/cuttingTriSurface.C
./surface/isoSurface/isoSurface.C
./surface/isoSurface/isoSurfaceCell.C
./surface/thresholdCellFaces/thresholdCellFaces.C
./surface/triSurfaceMesh/discreteSurface.C
./surfMeshSampler/surfMeshSampler/surfMeshSampler.C
./surfMeshSampler/surfMeshSamplers/surfMeshSamplers.C
./surfMeshSampler/plane/surfMeshPlaneSampler.C
./surfMeshSampler/triSurfaceMesh/surfMeshDiscreteSampler.C
./sampledSurface/sampledPatch/sampledPatch.C
./sampledSurface/sampledPatchInternalField/sampledPatchInternalField.C
./sampledSurface/sampledPlane/sampledPlane.C
./sampledSurface/isoSurface/sampledIsoSurface.C
./sampledSurface/isoSurface/sampledIsoSurfaceCell.C
./sampledSurface/distanceSurface/distanceSurface.C
./sampledSurface/sampledCuttingPlane/sampledCuttingPlane.C
./sampledSurface/sampledSurface/sampledSurface.C
./sampledSurface/sampledSurfaces/sampledSurfaces.C
./sampledSurface/sampledSurfaces/sampledSurfacesGrouping.C
./sampledSurface/sampledTriSurfaceMesh/sampledTriSurfaceMesh.C
./sampledSurface/sampledTriSurfaceMesh/sampledTriSurfaceMeshNormal.C
./sampledSurface/thresholdCellFaces/sampledThresholdCellFaces.C
./sampledSurface/sampledTriSurfaceToMesh/sampledTriSurfaceToMesh.C
./sampledSurface/writers/surfaceWriter.C
./sampledSurface/writers/ensight/ensightSurfaceWriter.C
./sampledSurface/writers/foam/foamSurfaceWriter.C
./sampledSurface/writers/nastran/nastranSurfaceWriter.C
./sampledSurface/writers/proxy/proxySurfaceWriter.C
./sampledSurface/writers/raw/rawSurfaceWriter.C
./sampledSurface/writers/starcd/starcdSurfaceWriter.C
./sampledSurface/writers/vtk/vtkSurfaceWriter.C
./sampledSurface/writers/vtkb/vtkbSurfaceWriter.C
./sampledSurface/writers/boundaryData/boundaryDataSurfaceWriter.C
./sampledSurface/readers/surfaceReader.C
./sampledSurface/readers/surfaceReaderNew.C
./sampledSurface/readers/ensight/ensightSurfaceReader.C
./graphField/writePatchGraph.C
./graphField/writeCellGraph.C
./graphField/makeGraph.C
./meshToMesh/meshToMesh.C
./meshToMesh/meshToMeshParallelOps.C
./meshToMesh/distributedWeightedFvPatchFieldMapper.C
./meshToMesh/calcMethod/meshToMeshMethod/meshToMeshMethod.C
./meshToMesh/calcMethod/meshToMeshMethod/meshToMeshMethodNew.C
./meshToMesh/calcMethod/cellVolumeWeight/cellVolumeWeightMethod.C
./meshToMesh/calcMethod/correctedCellVolumeWeight/correctedCellVolumeWeightMethod.C
./meshToMesh/calcMethod/direct/directMethod.C
./meshToMesh/calcMethod/mapNearest/mapNearestMethod.C
./meshToMesh0/meshToMesh0.C
./meshToMesh0/calculateMeshToMesh0Addressing.C
./meshToMesh0/calculateMeshToMesh0Weights.C
)

helyx_additional_includes(${TARGET_NAME}
    PUBLIC
        surfMesh-includes
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE finiteVolume
    PRIVATE meshTools
    PRIVATE conversion
    PRIVATE dynamicMesh
    PRIVATE lagrangian
    PRIVATE fileFormats
)

set_source_files_properties(
    ${CMAKE_CURRENT_LIST_DIR}/surface/cuttingTriSurface/cuttingTriSurface.C
    ${CMAKE_CURRENT_LIST_DIR}/sampledSet/array/arraySet.C
    ${CMAKE_CURRENT_LIST_DIR}/sampledSet/polyLine/polyLineSet.C
    ${CMAKE_CURRENT_LIST_DIR}/sampledSet/face/faceOnlySet.C
    ${CMAKE_CURRENT_LIST_DIR}/sampledSet/midPoint/midPointSet.C
    ${CMAKE_CURRENT_LIST_DIR}/sampledSet/midPointAndFace/midPointAndFaceSet.C

    # Redefinition of nearestEqOp, also defined in discreteSurface.C:65
    ${CMAKE_CURRENT_LIST_DIR}/sampledSurface/sampledTriSurfaceMesh/sampledTriSurfaceMesh.C

    PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON
)
