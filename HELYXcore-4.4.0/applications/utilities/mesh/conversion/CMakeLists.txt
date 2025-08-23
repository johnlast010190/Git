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

add_subdirectory(ansysToFoam)

add_subdirectory(cfx4ToFoam)

add_subdirectory(datToFoam)

add_subdirectory(fireToFoam)

add_subdirectory(fluent3DMeshToFoam)

add_subdirectory(fluentMeshToFoam)

add_subdirectory(foamMeshToFluent)

add_subdirectory(foamToFireMesh)

add_subdirectory(foamToStarMesh)

add_subdirectory(foamToSurface)

add_subdirectory(gambitToFoam)

add_subdirectory(gmshToFoam)

add_subdirectory(ideasUnvToFoam)

add_subdirectory(kivaToFoam)

add_subdirectory(mshToFoam)

add_subdirectory(netgenNeutralToFoam)

add_subdirectory(plot3dToFoam)

add_subdirectory(star4ToFoam)

add_subdirectory(tetgenToFoam)

add_subdirectory(vtkUnstructuredToFoam)

add_subdirectory(writeMeshObj)


