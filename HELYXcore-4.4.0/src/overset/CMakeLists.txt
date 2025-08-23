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

set(TARGET_NAME "overset")

add_helyx_library(${TARGET_NAME} SHARED
./cellCellStencil/cellCellStencil/cellCellStencil.C
./cellCellStencil/cellVolumeWeight/cellVolumeWeightCellCellStencil.C
./cellCellStencil/cellCellStencil/cellCellStencilObject.C
./cellCellStencil/inverseDistance/inverseDistanceCellCellStencil.C
./cellCellStencil/inverseDistance/waveMethod.C
./cellCellStencil/inverseDistance/meshToMeshData.C
./cellCellStencil/trackingInverseDistance/voxelMeshSearch.C
./cellCellStencil/trackingInverseDistance/trackingInverseDistanceCellCellStencil.C
./dynamicOversetFvMesh/dynamicOversetFvMesh.C
./fvMeshPrimitiveLduAddressing/fvMeshPrimitiveLduAddressing.C
./oversetPolyPatch/oversetPolyPatch.C
./oversetPolyPatch/oversetLduInterface.C
./oversetPolyPatch/oversetFvPatch.C
./oversetPolyPatch/oversetFvPatchFields.C
./oversetPolyPatch/oversetFvsPatchFields.C
./oversetPolyPatch/oversetGAMGInterface.C
./oversetPolyPatch/oversetGAMGInterfaceField.C
./oversetPolyPatch/oversetPointPatch.C
./oversetPolyPatch/oversetPointPatchFields.C
./oversetAdjustPhi/oversetAdjustPhi.C
./regionsToCell/regionsToCell.C
./lduPrimitiveProcessorInterface/lduPrimitiveProcessorInterface.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        fileFormats-includes
        surfMesh-includes
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE finiteVolume
    PRIVATE dynamicMesh
    PUBLIC dynamicFvMesh
    PRIVATE sampling
    PRIVATE meshTools
    PRIVATE decompositionMethods
    )
