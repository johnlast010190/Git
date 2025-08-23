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
    (c) 2016-2023 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "functionObjects/moleFractions/moleFractions.H"
#include "fluidMulticomponentThermo/fluidMulticomponentThermo.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(moleFractions, 0);
    addToRunTimeSelectionTable(functionObject, moleFractions, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::moleFractions::moleFractions
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    phaseName_(dict.lookupOrDefault<word>("phase", word::null)),
    X_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::moleFractions::~moleFractions()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::moleFractions::read(const dictionary& dict)
{
    return true;
}


bool Foam::functionObjects::moleFractions::execute()
{
    // Lookup or construct a multicomponent thermo. Lookup is likely to be
    // appropriate if this is being used at run-time. Construction if this is
    // being called as a one-off post-process.
    const word thermoName =
        IOobject::groupName(basicThermo::dictName, phaseName_);
    autoPtr<fluidMulticomponentThermo> thermoPtr
    (
        obr().foundObject<fluidMulticomponentThermo>(thermoName)
      ? autoPtr<fluidMulticomponentThermo>(nullptr)
      : fluidMulticomponentThermo::New(obr())
    );
    const fluidMulticomponentThermo& thermo =
        obr().lookupObject<fluidMulticomponentThermo>(thermoName);

    // Construct mole fraction fields corresponding to the mass fraction fields
    const PtrList<volScalarField>& Y = thermo.composition().Y();
    if (X_.empty())
    {
        X_.setSize(Y.size());

        forAll(Y, i)
        {
            X_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "X_" + Y[i].name(),
                        mesh_.time().timeName(),
                        obr()
                    ),
                    mesh_,
                    dimensionedScalar(dimless, 0)
                )
            );
        }
    }

    // Get the mixture molar mass
    const volScalarField W(thermo.W());

    // Calculate the mole fractions
    forAll(Y, i)
    {
        const dimensionedScalar Wi
        (
            "Wi",
            dimMass/dimMoles,
            thermo.composition().W(i)
        );

        X_[i] = Y[i]*W/Wi;
    }

    return true;
}


bool Foam::functionObjects::moleFractions::write()
{
    forAll(X_, i)
    {
        X_[i].write();
    }

    return true;
}


// ************************************************************************* //
