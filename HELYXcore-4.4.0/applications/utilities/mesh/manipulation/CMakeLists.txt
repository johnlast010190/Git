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

add_subdirectory(attachMesh)

add_subdirectory(autoPatch)

add_subdirectory(checkAndRemove)

add_subdirectory(checkMesh)

add_subdirectory(createBaffles)

add_subdirectory(createJump)

add_subdirectory(createNonConformalCouples)

add_subdirectory(createPatch)

add_subdirectory(deformedGeom)

add_subdirectory(flattenMesh)

add_subdirectory(helyxCheckMesh)

add_subdirectory(insideCells)

add_subdirectory(layerCalculator)

add_subdirectory(mergeMeshes)

add_subdirectory(mergeOrSplitBaffles)

add_subdirectory(meshOptimization)

add_subdirectory(mirrorMesh)

add_subdirectory(moveDynamicMesh)

add_subdirectory(moveEngineMesh)

add_subdirectory(moveMesh)

add_subdirectory(objToVTK)

add_subdirectory(orientFaceZone)

add_subdirectory(polyDualMesh)

add_subdirectory(refineMesh)

add_subdirectory(renumberMesh)

add_subdirectory(rotateMesh)

add_subdirectory(setSet)

add_subdirectory(setsToZones)

add_subdirectory(singleCellMesh)

add_subdirectory(splitMesh)

add_subdirectory(splitMeshRegions)

add_subdirectory(stitchMesh)

add_subdirectory(stretchMesh)

add_subdirectory(subsetMesh)

add_subdirectory(topoSet)

add_subdirectory(transformPoints)

add_subdirectory(zipUpMesh)


