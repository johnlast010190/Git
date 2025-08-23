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

set(TARGET_NAME "blockMesh")

add_helyx_library(${TARGET_NAME} SHARED
./blockVertices/blockVertex/blockVertex.C
./blockVertices/pointVertex/pointVertex.C
./blockVertices/projectVertex/projectVertex.C
./blockVertices/namedVertex/namedVertex.C
./blockEdges/bezier/bezier.C
./blockEdges/blockEdge/blockEdge.C
./blockEdges/lineDivide/lineDivide.C
./blockEdges/lineEdge/lineEdge.C
./blockEdges/polyLineEdge/polyLine.C
./blockEdges/polyLineEdge/polyLineEdge.C
./blockEdges/arcEdge/arcEdge.C
./blockEdges/BSplineEdge/BSpline.C
./blockEdges/BSplineEdge/BSplineEdge.C
./blockEdges/splineEdge/CatmullRomSpline.C
./blockEdges/splineEdge/splineEdge.C
./blockEdges/projectEdge/projectEdge.C
./blockEdges/projectCurveEdge/projectCurveEdge.C
./blockFaces/blockFace/blockFace.C
./blockFaces/projectFace/projectFace.C
./gradingDescriptor/gradingDescriptor.C
./gradingDescriptor/gradingDescriptors.C
./blockDescriptor/blockDescriptor.C
./blockDescriptor/blockDescriptorEdges.C
./blocks/block/block.C
./blocks/block/blockCreate.C
./blocks/namedBlock/namedBlock.C
./blockMesh/blockMesh.C
./blockMesh/blockMeshCreate.C
./blockMesh/blockMeshTopology.C
./blockMesh/blockMeshCheck.C
./blockMesh/blockMeshMerge.C
./blockMesh/blockMeshMergeFast.C
./blockMeshTools/blockMeshTools.C
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE meshTools
)
