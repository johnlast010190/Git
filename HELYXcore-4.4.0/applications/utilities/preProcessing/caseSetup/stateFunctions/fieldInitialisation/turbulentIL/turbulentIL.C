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

#include "turbulentIL.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fieldInitializations
{
    defineTypeNameAndDebug(turbulentIL, 0);
    addToRunTimeSelectionTable(fieldInit, turbulentIL, initMethod);

// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * //

scalar fieldInitializations::turbulentIL::calck(scalar I, scalar U)
{
    return 3.0/2.0 * Foam::sqr(U*I);
}

scalar turbulentIL::calcEpsilon(scalar I, scalar L, scalar U)
{
    scalar k = calck(I,U);
    return Foam::pow(0.09,(3.0/4.0)) * Foam::pow(k,(3.0/2.0)) / L;
}

scalar turbulentIL::calcNut(scalar I, scalar L, scalar U)
{
    scalar k = calck(I,U);
    return Foam::pow(0.09,(1.0/4.0))*Foam::sqrt(k)*L;
}

scalar turbulentIL::calcOmega(scalar I, scalar L, scalar U)
{
    scalar k = calck(I,U);
    return Foam::pow(0.09,(-1.0/4.0)) * Foam::sqrt(k) / L;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentIL::turbulentIL
(
    const fvMesh& mesh,
    const fvSolutionRegistry& localDb,
    const dictionary& fieldDict,
    const word& fN
)
:
    fieldInit(mesh, localDb, fieldDict, fN)
{
    // dependent on velcotiy field to be initialised first
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void turbulentIL::correct()
{

    //if already initialised stop
    if (initialised())
    {
        return;
    }

    fieldInit::initMsg();

    //always scalar field
    volScalarField& f
    (
        const_cast<volScalarField&>
        (localDb().lookupObject<volScalarField>(name()))
    );

    scalar I = initDict().lookup<scalar>("I");
    scalar L = initDict().lookup<scalar>("L");
    scalar Uref = initDict().lookup<scalar>("Uref");

    word group = IOobject::group(name());

    word UName = word
    (
        initDict().lookupOrDefault<word>("U", IOobject::groupName("U", group))
    );

    const volVectorField* UPtr(nullptr);

    if (localDb().foundObject<volVectorField>(UName))
    {
        UPtr = &(localDb().lookupObject<volVectorField>(UName));
    }
    else
    {
        FatalErrorInFunction
            << "cannot find velocity field " << UName
            << " to initialize turbulent values"
            << exit(FatalError);
    }

    const volVectorField& U = *UPtr;

    if (name() == "k" || name() == "kl")
    {
        forAll(f, fI)
        {
            scalar Uturb = Uref;
            if (mag(U[fI]) > SMALL)
            {
                Uturb = mag(U[fI]);
            }
            f[fI] = calck(I, Uturb);
        }
    }
    else if (name() == "epsilon")
    {
        forAll(f, fI)
        {
            scalar Uturb = Uref;
            if (mag(U[fI]) > SMALL)
            {
                Uturb = mag(U[fI]);
            }
            f[fI] = calcEpsilon(I, L, Uturb);
        }
    }
    else if (name() == "nuTilda")
    {
        forAll(f, fI)
        {
            scalar Uturb = Uref;
            if (mag(U[fI]) > SMALL)
            {
                Uturb = mag(U[fI]);
            }
            f[fI] = calcNut(I, L, Uturb);
        }
    }
    else if (name() == "nut")
    {
        forAll(f, fI)
        {
            scalar Uturb = Uref;
            if (mag(U[fI]) > SMALL)
            {
                Uturb = mag(U[fI]);
            }
            f[fI] = calcNut(I, L, Uturb);
        }
    }
    else if (name() == "mut")
    {
        scalar rho = initDict().lookup<scalar>("rho");

        forAll(f, fI)
        {
            scalar Uturb = Uref;

            if (mag(U[fI]) > SMALL)
            {
                Uturb = mag(U[fI]);
            }

            f[fI] = rho*calcNut(I, L, Uturb);
        }
    }
    else if (name() == "omega")
    {
        forAll(f, fI)
        {
            scalar Uturb = Uref;

            if (mag(U[fI]) > SMALL)
            {
                Uturb = mag(U[fI]);
            }

            f[fI] = calcOmega(I, L, Uturb);
        }
    }

    // the field has been initialised
    initialised() = true;

    initCoupledBoundaries();

    initTurbBoundaries();
}

} //end namespace fieldInitialisations
} //end namespace FOAM


// ************************************************************************* //
