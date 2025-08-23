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
    (c) 2011-2023 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "meshes/pointMesh/pointMesh.H"
#include "meshes/polyMesh/globalMeshData/globalMeshData.H"
#include "fields/GeometricFields/pointFields/pointFields.H"
#include "meshes/pointMesh/pointPatches/facePointPatch/facePointPatch.H"
#include "fields/GeometricFields/GeometricField/MapGeometricFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointMesh, 0);
}

const Foam::HashSet<Foam::word> Foam::pointMesh::geometryFields;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointMesh::pointMesh(const polyMesh& pMesh)
:
    MeshObject<polyMesh, Foam::PatchMeshObject, pointMesh>(pMesh),
    GeoMesh<polyMesh>(pMesh),
    boundary_(*this, pMesh.boundaryMesh())
{
    if (debug)
    {
        Pout<< "pointMesh::pointMesh(const polyMesh&): "
            << "Constructing from polyMesh " << pMesh.name()
            << endl;
    }

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointMesh::~pointMesh()
{
    if (debug)
    {
        Pout<< "~pointMesh::pointMesh()"
            << endl;
        error::printStack(Pout);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::pointMesh::movePoints()
{
    if (debug)
    {
        Pout<< "pointMesh::movePoints(const pointField&): "
            << "Moving points." << endl;
    }

    boundary_.movePoints(GeoMesh<polyMesh>::mesh_.points());

    return true;
}


void Foam::pointMesh::topoChange(const polyTopoChangeMap& map)
{
    if (debug)
    {
        Pout<< "pointMesh::topoChange(const polyTopoChangeMap&): "
            << "Topology change." << endl;
        Pout<< endl;
    }
    boundary_.topoChange();
}


void Foam::pointMesh::mapMesh(const polyMeshMap& map)
{
    if (debug)
    {
        Pout<< "pointMesh::mapMesh(const polyMeshMap&): "
            << "Mesh mapping." << endl;
        Pout<< endl;
    }
    boundary_.topoChange();
}


void Foam::pointMesh::distribute(const polyDistributionMap& map)
{
    if (debug)
    {
        Pout<< "pointMesh::distribute(const polyDistributionMap&): "
            << "Distribute." << endl;
        Pout<< endl;
    }
    boundary_.topoChange();
}


void Foam::pointMesh::reorderPatches
(
    const labelUList& newToOld,
    const bool validBoundary
)
{
    if (debug)
    {
        Pout<< "pointMesh::reorderPatches( const labelUList&, const bool): "
            << "Updating for reordered patches." << endl;
        Pout<< endl;
    }

    boundary_.shuffle(newToOld, validBoundary);

    objectRegistry& db = const_cast<objectRegistry&>(thisDb());
    ReorderPatchFields<pointScalarField>(db, newToOld);
    ReorderPatchFields<pointVectorField>(db, newToOld);
    ReorderPatchFields<pointSphericalTensorField>(db, newToOld);
    ReorderPatchFields<pointSymmTensorField>(db, newToOld);
    ReorderPatchFields<pointTensorField>(db, newToOld);
}


void Foam::pointMesh::addPatch(const label patchi)
{
    if (debug)
    {
        Pout<< "pointMesh::addPatch(const label): "
            << "Adding patch at " << patchi << endl;
        Pout<< endl;
    }

    const polyBoundaryMesh& pbm = mesh().boundaryMesh();
    if (pbm.size() != boundary_.size())
    {
        FatalErrorInFunction << "Problem :"
            << " pointBoundaryMesh size :" << boundary_.size()
            << " polyBoundaryMesh size :" << pbm.size()
            << exit(FatalError);
    }

    boundary_.set(patchi, facePointPatch::New(pbm[patchi], boundary_).ptr());

    objectRegistry& db = const_cast<objectRegistry&>(thisDb());
    const dictionary d;
    const word patchFieldType("calculated");

    AddPatchFields<pointScalarField>(db, patchi, d, patchFieldType, Zero);
    AddPatchFields<pointVectorField>(db, patchi, d, patchFieldType, Zero);
    AddPatchFields<pointSphericalTensorField>
    (
        db,
        patchi,
        d,
        patchFieldType,
        Zero
    );
    AddPatchFields<pointSymmTensorField>(db, patchi, d, patchFieldType, Zero);
    AddPatchFields<pointTensorField>(db, patchi, d, patchFieldType, Zero);
}


void Foam::pointMesh::reset()
{
    if (debug)
    {
        Pout<< "pointMesh::reset(): "
            << "Mesh reset." << endl;
        Pout<< endl;
    }
    boundary_.reset();
}


// ************************************************************************* //
