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

add_subdirectory(surfaceAdd)

add_subdirectory(surfaceCheck)

add_subdirectory(surfaceClean)

add_subdirectory(surfaceCoarsen)

add_subdirectory(surfaceConvert)

add_subdirectory(surfaceFeatureConvert)

add_subdirectory(surfaceFeatureExtract)

add_subdirectory(surfaceFind)

add_subdirectory(surfaceFix)

add_subdirectory(surfaceHookUp)

add_subdirectory(surfaceInertia)

add_subdirectory(surfaceInflate)

add_subdirectory(surfaceInsidePoint)

add_subdirectory(surfaceIntersection)

add_subdirectory(surfaceLambdaMuSmooth)

add_subdirectory(surfaceMeshConvert)

add_subdirectory(surfaceMeshConvertTesting)

add_subdirectory(surfaceMeshExport)

add_subdirectory(surfaceMeshImport)

add_subdirectory(surfaceMeshInfo)

add_subdirectory(surfaceMeshTriangulate)

add_subdirectory(surfaceOrient)

add_subdirectory(surfacePatch)

add_subdirectory(surfacePointMerge)

add_subdirectory(surfaceRedistributePar)

add_subdirectory(surfaceRefineRedGreen)

add_subdirectory(surfaceRotationAxis)

add_subdirectory(surfaceSplitByPatch)

add_subdirectory(surfaceSplitByTopology)

add_subdirectory(surfaceSplitNonManifolds)

add_subdirectory(surfaceSubset)

add_subdirectory(surfaceToPatch)

add_subdirectory(surfaceTransformPoints)

if (${OPENVDB_FOUND} AND ${OPENCASCADE_FOUND} AND ("${HELYX_LABEL_SIZE}" STREQUAL "32"))
    add_subdirectory(helyxVDBSurfaceWrap)
else ()
    add_dummy_target(helyxVDBSurfaceWrap OpenCASCADE)
endif ()

