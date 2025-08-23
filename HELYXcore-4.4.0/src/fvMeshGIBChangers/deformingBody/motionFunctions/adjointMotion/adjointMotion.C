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
    (c) 2011 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "deformingBody/motionFunctions/adjointMotion/adjointMotion.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshGIBChangers
{
    defineTypeNameAndDebug(adjointMotion, 0);

    addToRunTimeSelectionTable
    (
        motionFunction,
        adjointMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::fvMeshGIBChangers::adjointMotion::normalVelocity
(
    const label& pI
) const
{
    tmp<scalarField > Unt(new scalarField(mesh_.boundary()[pI].size()));
    scalarField& Un = Unt.ref();

    volScalarField& G =
        const_cast<volScalarField&>(mesh_.lookupObject<volScalarField>("G"));

    forAll(Un, fI)
    {
        Un[fI] = G.boundaryField()[pI][fI];
        Un[fI] = min(Un[fI],0.3);
        Un[fI] = max(Un[fI],-0.3);
    }

    return Unt;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshGIBChangers::adjointMotion::adjointMotion
(
    const fvMesh& mesh,
    const dictionary& DBMFCoeffs,
    const Time& runTime
)
:
    motionFunction(mesh, DBMFCoeffs, runTime)
{
    read(DBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::fvMeshGIBChangers::adjointMotion::~adjointMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField>
Foam::fvMeshGIBChangers::adjointMotion::boundaryVelocity
(
    const label& pMaster
) const
{
    tmp<vectorField > tVelocityBoundary
    (
        new vectorField(mesh_.boundary()[pMaster].size(), vector::zero)
    );

    return tVelocityBoundary;
}


Foam::tmp<Foam::vectorField>
Foam::fvMeshGIBChangers::adjointMotion::interfaceVelocity
(
    const label& pMaster
) const
{
    tmp<scalarField> Unt = normalVelocity(pMaster);

    tmp<vectorField> tSpeedInter(new vectorField(Unt->size(), vector::zero));
    vectorField& speedInter = tSpeedInter.ref();

    speedInter = Unt()*mesh_.boundary()[pMaster].nf()();

    return tSpeedInter;
}


// ************************************************************************* //
