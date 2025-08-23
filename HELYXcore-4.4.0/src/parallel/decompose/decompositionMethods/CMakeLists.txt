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

set(TARGET_NAME "decompositionMethods")

add_helyx_library(${TARGET_NAME} SHARED
./decompositionMethod/decompositionMethod.C
./geomDecomp/geomDecomp.C
./simpleGeomDecomp/simpleGeomDecomp.C
./hierarchGeomDecomp/hierarchGeomDecomp.C
./manualDecomp/manualDecomp.C
./metisLikeDecomp/metisLikeDecomp.C
./multiLevelDecomp/multiLevelDecomp.C
./structuredDecomp/structuredDecomp.C
./noDecomp/noDecomp.C
./decompositionConstraints/decompositionConstraint/decompositionConstraint.C
./decompositionConstraints/preserveBaffles/preserveBafflesConstraint.C
./decompositionConstraints/preserveFaceZones/preserveFaceZonesConstraint.C
./decompositionConstraints/preservePatches/preservePatchesConstraint.C
./decompositionConstraints/singleProcessorFaceSets/singleProcessorFaceSetsConstraint.C
./decompositionConstraints/refinementHistory/refinementHistoryConstraint.C
./decompositionConstraints/singleProcessorPatches/singleProcessorPatchesConstraint.C
)

# Oh god please introduce an include directory if your'e going to put this stuff in a shared library.
target_include_directories(${TARGET_NAME}
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/decompositionMethod>
)


helyx_link_libraries(${TARGET_NAME}
    PRIVATE meshTools
    PRIVATE dynamicMesh
    PRIVATE finiteVolume
)
