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

set(TARGET_NAME "conversion")

add_helyx_library(${TARGET_NAME} SHARED
./common/reader/meshReader.C
./common/reader/meshReaderAux.C
./common/reader/calcPointCells.C
./common/reader/createPolyCells.C
./common/reader/createPolyBoundary.C
./common/writer/meshWriter.C
./common/tables/boundaryRegion.C
./common/tables/cellTable.C
./ensight/mesh/ensightMesh.C
./ensight/mesh/ensightMeshIO.C
./ensight/mesh/ensightMeshOptions.C
./ensight/part/ensightPart.C
./ensight/part/ensightPartCells.C
./ensight/part/ensightPartFaces.C
./ensight/part/ensightParts.C
./fire/FIREMeshReader.C
./fire/FIREMeshWriter.C
./fire/checkFireEdges.C
./mergedCyclic/mergedCyclicPolyPatch.C
./mergedCyclic/polyMeshUnMergeCyclics.C
./starcd/STARCDMeshReader.C
./starcd/STARCDMeshWriter.C
./polyDualMesh/polyDualMesh.C
./vtk/output/foamVtkInternalWriter.C
./vtk/output/foamVtkPatchWriter.C
./vtk/output/foamVtkSurfaceFieldWriter.C
)

helyx_additional_includes(${TARGET_NAME}
    PUBLIC
        fileFormats-includes
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE finiteVolume
    PRIVATE finiteArea
    PRIVATE meshTools
    PRIVATE fileFormats
)
