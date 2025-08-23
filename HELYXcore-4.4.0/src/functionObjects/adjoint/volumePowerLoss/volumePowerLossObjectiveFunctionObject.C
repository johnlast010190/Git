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
    (c) 2019 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "volumePowerLoss/volumePowerLossObjectiveFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(volumePowerLossObjectiveFunctionObject, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        volumePowerLossObjectiveFunctionObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::volumePowerLossObjectiveFunctionObject::
volumePowerLossObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict,
    bool useAdjointFileFormat
)
:
    objectiveFunctionObject(runTime, const_cast<dictionary&>(objectiveDict)),
    AreaObjPatchInlet_(objectivePatchInletArea()),
    AreaObjPatchOutlet_(objectivePatchOutletArea()),
    adiabatic_(objectiveDict.lookupOrDefault<Switch>("adiabatic", false)),
    stressWork_(objectiveDict.lookupOrDefault<Switch>("stressWork", false)),
    printLossesSeparately_
    (
        objectiveDict.lookupOrDefault<Switch>("printLossesSeparately", false)
    )
{
    createFiles(useAdjointFileFormat);
}


Foam::functionObjects::volumePowerLossObjectiveFunctionObject::
volumePowerLossObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict
)
:
    volumePowerLossObjectiveFunctionObject(name, runTime, objectiveDict, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::volumePowerLossObjectiveFunctionObject::
~volumePowerLossObjectiveFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::functionObjects::volumePowerLossObjectiveFunctionObject::read
(
    const dictionary& dict
)
{
    objectiveFunctionObject::read(dict);

    AreaObjPatchInlet_ = objectivePatchInletArea();
    AreaObjPatchOutlet_ = objectivePatchOutletArea();
    adiabatic_ = dict.lookupOrDefault<Switch>("adiabatic", false);
    stressWork_ = dict.lookupOrDefault<Switch>("stressWork", false);

    return true;
}


bool
Foam::functionObjects::volumePowerLossObjectiveFunctionObject::execute()
{
    scalar powerLoss = 0;
    scalar viscousLosses = 0;
    scalar porosityLosses = 0;

    tmp<volSymmTensorField> tdevRhoReff = devRhoReff();

    viscousLosses =
        fvc::domainIntegrate
        (
          - (tdevRhoReff && fvc::grad(U()))
        ).value();

    bool isUsingPorosity =
        mesh_.foundObject<volScalarField>("porosity");

    if (isUsingPorosity)
    {
        const volScalarField& porosity =
            mesh_.lookupObject<volScalarField>("porosity");

        porosityLosses =
            fvc::domainIntegrate
            (
                rho()*porosity*magSqr(U())
            ).value();
    }

    powerLoss = viscousLosses + porosityLosses;
    objectiveValue_ = powerLoss;

    Info<< type() << " " << name() << " execute:" << nl
        << "Power losses = " << objectiveValue_ << " [W]" << nl;

    if (printLossesSeparately_)
    {
        Info<< "Viscous Losses = " << viscousLosses << " [W], "
            << "Porosity Losses = " << porosityLosses << " [W]";
    }

    Info<< endl;

    return true;
}


bool
Foam::functionObjects::volumePowerLossObjectiveFunctionObject::write()
{
    writeOutputToFile();

    return true;
}


// ************************************************************************* //
