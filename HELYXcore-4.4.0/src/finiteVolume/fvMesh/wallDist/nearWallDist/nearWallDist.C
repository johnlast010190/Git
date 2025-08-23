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
    (c) 2017-2023 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMesh/wallDist/nearWallDist/nearWallDist.H"
#include "fvMesh/fvMesh.H"
#include "patchDist/patchDistFuncs/patchDistFuncs.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "fvMesh/fvPatches/derived/mapped/mappedWallFvPatch.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nearWallDist, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::nearWallDist::resize()
{
    y_.setSize(mesh().boundary().size());

    forAll(y_, patchi)
    {
        y_.set
        (
            patchi,
            fvPatchField<scalar>::New
            (
                calculatedFvPatchScalarField::typeName,
                mesh().boundary()[patchi],
                volScalarField::Internal::null()
            )
        );
    }
}


void Foam::nearWallDist::correct()
{
    if (mesh().topoChanged())
    {

    }

    patchDistFuncs::correctBoundaryFaceFaceCells
    (
        mesh(),
        mesh().boundaryMesh().findPatchIDs<wall>(),
        y_
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nearWallDist::nearWallDist(const Foam::fvMesh& mesh)
:
    MeshObject<fvMesh, Foam::UpdateableMeshObject, nearWallDist>(mesh),
    y_
    (
        mesh.boundary(),
        volScalarField::Internal::null(),
        calculatedFvPatchScalarField::typeName
    )
{
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nearWallDist::~nearWallDist()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::nearWallDist::movePoints()
{
    resize();
    correct();
    return true;
}


void Foam::nearWallDist::topoChange(const polyTopoChangeMap& map)
{
    resize();
    correct();
}


void Foam::nearWallDist::mapMesh(const polyMeshMap& map)
{
    resize();
    correct();
}


void Foam::nearWallDist::distribute(const polyDistributionMap& map)
{
    resize();
    correct();
}


// ************************************************************************* //
