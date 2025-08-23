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
    (c) 2016 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "relativeHumidity.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "basicThermo/basicThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fieldInitializations
{
    defineTypeNameAndDebug(relativeHumidity, 0);
    addToRunTimeSelectionTable(fieldInit, relativeHumidity, initMethod);
}
}

// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldInitializations::relativeHumidity::relativeHumidity
(
    const fvMesh& mesh,
    const fvSolutionRegistry& localDb,
    const dictionary& fieldDict,
    const word& fN
)
:
    fieldInit(mesh, localDb, fieldDict, fN)
{
    //requires temperature field
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fieldInitializations::relativeHumidity::correct()
{
    //if already initialized stop
    if (initialised())
    {
        return;
    }

    fieldInit::initMsg();

    word TName = initDict().lookupOrDefault<word>("T", "T");
    tmp<volScalarField> TAbs = basicThermo::TAbsIfFound(localDb(), TName);

    scalar pAbs = initDict().lookupOrDefault<scalar>("Pabs", 1e5);
    scalar Mvap = initDict().lookupOrDefault<scalar>("Mvap", 18.02);
    scalar Mmix = initDict().lookupOrDefault<scalar>("Mmix", 28.96);
    scalar relativeHumidity
        = initDict().lookup<scalar>("relativeHumidity");

    //hardcoded saturation pressure for water vapour

    scalarField Ppartial
    (
        (2337 * Foam::exp(6879*(-1/TAbs().primitiveField() + 1/293.15)
        - 5.031*Foam::log(TAbs().primitiveField()/293.15)))
        * relativeHumidity
    );

    volScalarField& f
        = const_cast<volScalarField&>
        (localDb().lookupObject<volScalarField>(name()));

    f.primitiveFieldRef()
        = Ppartial * Mvap / (Ppartial * Mvap + (pAbs - Ppartial)*Mmix);

    // the field has been initialized
    initialised() = true;

    initCoupledBoundaries();
}



// ************************************************************************* //
