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

set(TARGET_NAME "surfMesh")

add_helyx_library(${TARGET_NAME} SHARED
./surfZone/surfZone/surfZone.C
./surfZone/surfZone/surfZoneIOList.C
./MeshedSurfaceAllocator/MeshedSurfaceIOAllocator.C
./MeshedSurface/MeshedSurfaceCore.C
./MeshedSurface/MeshedSurfaces.C
./UnsortedMeshedSurface/UnsortedMeshedSurfaces.C
./MeshedSurfaceProxy/MeshedSurfaceProxys.C
./mergedSurf/mergedSurf.C
./polySurface/polySurface.C
./polySurface/polySurfaceClear.C
./polySurface/polySurfaceIO.C
./polySurface/fields/polySurfaceFields.C
./polySurface/fields/polySurfacePointFields.C
./surfaceRegistry/surfaceRegistry.C
./meshToSurface/meshToSurface.C
./surfMesh/surfMesh.C
./surfMesh/surfMeshClear.C
./surfMesh/surfMeshIO.C
./surfFields/surfFields/surfFields.C
./surfFields/surfPointFields/surfPointFields.C
./surfaceFormats/surfaceFormatsCore.C
./surfaceFormats/ac3d/AC3DsurfaceFormatCore.C
./surfaceFormats/ac3d/AC3DsurfaceFormatRunTime.C
./surfaceFormats/fire/FLMAsurfaceFormatRunTime.C
./surfaceFormats/gts/GTSsurfaceFormatRunTime.C
./surfaceFormats/nas/NASsurfaceFormatRunTime.C
./surfaceFormats/obj/OBJsurfaceFormatRunTime.C
./surfaceFormats/obj/OBJstream.C
./surfaceFormats/off/OFFsurfaceFormatRunTime.C
./surfaceFormats/smesh/SMESHsurfaceFormatRunTime.C
./surfaceFormats/starcd/STARCDsurfaceFormatCore.C
./surfaceFormats/starcd/STARCDsurfaceFormatRunTime.C
./surfaceFormats/stl/STLsurfaceFormatRunTime.C
./surfaceFormats/tri/TRIsurfaceFormatCore.C
./surfaceFormats/tri/TRIsurfaceFormatRunTime.C
./surfaceFormats/vtk/VTKsurfaceFormatCore.C
./surfaceFormats/vtk/VTKsurfaceFormatRunTime.C
./surfaceFormats/vtp/VTPsurfaceFormatCore.C
./surfaceFormats/vtp/VTPsurfaceFormatRunTime.C
./surfaceFormats/x3d/X3DsurfaceFormatCore.C
./surfaceFormats/x3d/X3DsurfaceFormatRunTime.C
./triSurface/triSurface.C
./triSurface/triSurfaceAddressing.C
./triSurface/stitchTriangles.C
./triSurface/fields/triSurfaceFields.C
./triSurface/patches/geometricSurfacePatch.C
./triSurface/patches/surfacePatch.C
./triSurface/interfaces/STL/writeSTL.C
./triSurface/interfaces/STL/readSTL.C
./triSurface/interfaces/GTS/writeGTS.C
./triSurface/interfaces/GTS/readGTS.C
./triSurface/interfaces/OBJ/readOBJ.C
./triSurface/interfaces/OBJ/writeOBJ.C
./triSurface/interfaces/SMESH/writeSMESH.C
./triSurface/interfaces/OFF/readOFF.C
./triSurface/interfaces/OFF/writeOFF.C
./triSurface/interfaces/TRI/writeTRI.C
./triSurface/interfaces/TRI/readTRI.C
./triSurface/interfaces/AC3D/readAC.C
./triSurface/interfaces/AC3D/writeAC.C
./triSurface/interfaces/VTK/readVTK.C
./triSurface/interfaces/VTK/writeVTK.C
./triSurface/interfaces/NAS/readNAS.C
./triSurface/interfaces/NAS/writeNAS.C
./triSurface/interfaces/INP/readINP.C
./triSurface/interfaces/INP/writeINP.C
./triSurface/interfaces/NTL/readNTL.C
./triSurface/interfaces/BTS/readBTS.C
)

if("${OPENCASCADE_FOUND}")
    target_sources(${TARGET_NAME}
        PRIVATE ./triSurface/interfaces/CAD/readCAD.C
    )
    helyx_link_libraries(${TARGET_NAME}
        PRIVATE ${THIRDPARTY_OPENCASCADE}
    )
    set_source_files_properties(
        ${CMAKE_CURRENT_LIST_DIR}/triSurface/interfaces/CAD/readCAD.C
        PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON
    )
endif()

if("${ALEMBIC_FOUND}")
    target_sources(${TARGET_NAME}
        PRIVATE ./triSurface/interfaces/ALEMBIC/readALEMBIC.C
    )
    helyx_link_libraries(${TARGET_NAME}
        PRIVATE ${THIRDPARTY_ALEMBIC}
    )
    set_source_files_properties(
        ${CMAKE_CURRENT_LIST_DIR}/triSurface/interfaces/ALEMBIC/readALEMBIC.C
        PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON
    )
else()
    target_sources(${TARGET_NAME}
        PRIVATE ./triSurface/interfaces/ALEMBIC/dummyReadALEMBIC.C
    )
endif()

helyx_link_libraries(${TARGET_NAME}
    PRIVATE fileFormats
)

set_source_files_properties(
    ${CMAKE_CURRENT_LIST_DIR}/triSurface/interfaces/VTK/writeVTK.C
    ${CMAKE_CURRENT_LIST_DIR}/surfaceFormats/vtp/VTPsurfaceFormatRunTime.C
    PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON
)
