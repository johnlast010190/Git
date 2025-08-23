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
    (c) 2018-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "deformingBody/motionFunctions/staticMotion/staticMotion.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshGIBChangers
{
    defineTypeNameAndDebug(staticMotion, 0);
    addToRunTimeSelectionTable
    (
        motionFunction,
        staticMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::fvMeshGIBChangers::staticMotion::normalVelocity
(
    const label& pI
) const
{
    tmp<scalarField> Unt(new scalarField(mesh_.boundary()[pI].size(), 0));

    return Unt;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshGIBChangers::staticMotion::staticMotion
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

Foam::fvMeshGIBChangers::staticMotion::~staticMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField>
Foam::fvMeshGIBChangers::staticMotion::boundaryVelocity
(
    const label& pMaster
) const
{
    tmp<vectorField> tVelocityBoundary
    (
        new vectorField(mesh_.boundary()[pMaster].size(), vector::zero)
    );

    return tVelocityBoundary;
}


Foam::tmp<Foam::vectorField>
Foam::fvMeshGIBChangers::staticMotion::interfaceVelocity
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
