/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
|    o     o     |  ENGYS Ltd. <http://engys.com/>
|       o        |
\*---------------------------------------------------------------------------
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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "surfMesh/surfMesh.H"
#include "meshes/polyMesh/globalMeshData/globalMeshData.H"
#include "include/demandDrivenData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::surfMesh::removeZones()
{
    if (debug)
    {
        InfoInFunction << "Removing surface zones." << endl;
    }

    // Remove the surface zones
    storedZones().clear();

    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfMesh::clearGeom()
{
    if (debug)
    {
        InfoInFunction << "Clearing geometric data" << endl;
    }

    MeshReference::clearGeom();
}


void Foam::surfMesh::clearAddressing()
{
    if (debug)
    {
        InfoInFunction << "clearing topology" << endl;
    }

    MeshReference::clearPatchMeshAddr();
}


void Foam::surfMesh::clearOut()
{
    MeshReference::clearOut();

    clearGeom();
    clearAddressing();
}


// ************************************************************************* //
