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
    (c) 2016-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "buoyant.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fieldInitializations
{
    defineTypeNameAndDebug(buoyant, 0);
    addToRunTimeSelectionTable(fieldInit, buoyant, initMethod);
}
}

// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldInitializations::buoyant::buoyant
(
    const fvMesh& mesh,
    const fvSolutionRegistry& localDb,
    const dictionary& fieldDict,
    const word& fN
)
:
    fieldInit(mesh, localDb, fieldDict, fN)
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fieldInitializations::buoyant::correct()
{
    //if already initialised stop
    if (initialised())
    {
        return;
    }

    fieldInit::initMsg();

    //valid for p only
    volScalarField& f
    (
        const_cast<volScalarField&>
        (localDb().lookupObject<volScalarField>(name()))
    );

    const label refCell = initDict().lookup<label>("refCell");
    const scalar refValue = initDict().lookup<scalar>("refValue");

    scalar rho = initDict().lookup<scalar>("rho");

    uniformDimensionedVectorField g
    (
        IOobject
        (
            "g",
            mesh().time().constant(),
            localDb(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    scalar refPoint = (g.value() & mesh().C()[refCell]);

    if (!Pstream::master())
    {
        refPoint = -GREAT;
    }

    reduce(refPoint, maxOp<scalar>());

    if (f.dimensions() == sqr(dimLength)/sqr(dimTime))
    {
        f.primitiveFieldRef() = (refValue * rho
            + rho*((g.value() & mesh().C().primitiveField()) - refPoint))
            /rho;
    }
    else
    {
        f.primitiveFieldRef() = refValue
            + rho*((g.value() & mesh().C().primitiveField()) - refPoint);
    }

    // the field has been initialised
    initialised() = true;

    initCoupledBoundaries();
}



// ************************************************************************* //
