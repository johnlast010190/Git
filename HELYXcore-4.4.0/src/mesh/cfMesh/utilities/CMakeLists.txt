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

add_subdirectory(FLMAToSurface)

add_subdirectory(FMSToSurface)

add_subdirectory(FMSToVTK)

add_subdirectory(checkSurfaceMesh)

add_subdirectory(copySurfaceParts)

add_subdirectory(extrudeEdgesInto2DSurface)

add_subdirectory(generateBoundaryLayers)

add_subdirectory(importSurfaceAsSubset)

add_subdirectory(improveMeshQuality)

add_subdirectory(improveSymmetryPlanes)

add_subdirectory(mergeSurfacePatches)

add_subdirectory(meshToFPMA)

add_subdirectory(patchesToSubsets)

add_subdirectory(preparePar)

add_subdirectory(removeSurfaceFacets)

add_subdirectory(scaleMesh)

add_subdirectory(scaleSurfaceMesh)

add_subdirectory(subsetToPatch)

add_subdirectory(surfaceFeatureEdges)

add_subdirectory(surfaceGenerateBoundingBox)

add_subdirectory(surfaceToFMS)


