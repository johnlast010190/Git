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
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "equationOfState.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "cfdTools/general/include/fvCFD.H"
#include "rhoThermo/rhoThermo.H"
#include "stateFunction/stateFunction.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fieldInitializations
{
    defineTypeNameAndDebug(equationOfState, 0);
    addToRunTimeSelectionTable(fieldInit, equationOfState, initMethod);
}
}

// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldInitializations::equationOfState::equationOfState
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

void Foam::fieldInitializations::equationOfState::correct()
{
    //if already initialised stop
    if (initialised())
    {
        return;
    }

    fieldInit::initMsg();

    // Calculate and write density

    // Protect T fields from modification by thermo
    PtrList<volScalarField> savedFields;
    HashTable<const volScalarField*> fields
    (
        localDb().lookupClass<volScalarField>()
    );
    forAllConstIters(fields, fieldIter)
    {
        if (IOobject::member(fieldIter.key()) == "T")
        {
            savedFields.append(new volScalarField(*fieldIter()));
        }
    }

    Info<< "Reading thermophysical properties\n" << endl;
    basicThermo* thermoPtr = &basicThermo::lookupOrCreate(localDb());

    volScalarField& rho = localDb().lookupObjectRef<volScalarField>(name());
    rho.forceAssign(thermoPtr->rho());

    // Restore original T fields
    forAll(savedFields, fieldi)
    {
        const_cast<volScalarField *>
        (
            fields[savedFields[fieldi].name()]
        )->forceAssign
        (
            savedFields[fieldi]
        );
    }

    // the field has been initialised
    initialised() = true;
}



// ************************************************************************* //
