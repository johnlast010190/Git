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
    (c) 2023 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "functionObjects/massFractions/massFractions.H"
#include "fluidMulticomponentThermo/fluidMulticomponentThermo.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(massFractions, 0);
    addToRunTimeSelectionTable(functionObject, massFractions, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::massFractions::massFractions
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    phaseName_(dict.lookupOrDefault<word>("phase", word::null)),
    Y_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::massFractions::~massFractions()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::massFractions::read(const dictionary& dict)
{
    return true;
}


bool Foam::functionObjects::massFractions::execute()
{
    // A multicomponent thermo should not exist for this context, as the
    // following would create a different, inconsistent definition of the
    // composition of an existing thermo model
    const word thermoName =
        IOobject::groupName(basicThermo::dictName, phaseName_);
    if (obr().foundObject<fluidMulticomponentThermo>(thermoName))
    {
        FatalErrorInFunction
            << "Cannot create mass fractions. Mass fractions already exist "
            << "within the multicomponent thermo model, \"" << thermoName
            << "\"." << exit(FatalError);
    }

    // Read Ydefault if it exists
    IOobject YdefaultIo
    (
        "Ydefault",
        time_.timeName(),
        mesh_,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );
    const bool YdefaultIoExists = YdefaultIo.headerOk();
    const volScalarField Ydefault
    (
        YdefaultIo,
        mesh_,
        dimensionedScalar(dimless, 0)
    );

    // Back up Ydefault if it exists
    if (YdefaultIoExists)
    {
        fileHandler().cp
        (
            YdefaultIo.filePath(),
            YdefaultIo.path(time_.instance())
           /"molesToMassFractions-" + YdefaultIo.name()
        );
    }

    // Write Ydefault out again, but limited to a value of small. This prevents
    // errors in construction of the thermo if no non-default fractions exist
    {
        volScalarField YdefaultLim(Ydefault);
        YdefaultLim.max(SMALL);
        YdefaultLim.write();
    }

    // Construct a multicomponent thermo
    autoPtr<fluidMulticomponentThermo> thermoPtr
    (
        obr().foundObject<fluidMulticomponentThermo>(thermoName)
      ? autoPtr<fluidMulticomponentThermo>(nullptr)
      : fluidMulticomponentThermo::New(obr())
    );
    fluidMulticomponentThermo& thermo =
        obr().lookupObjectRef<fluidMulticomponentThermo>(thermoName);
    PtrList<volScalarField>& Y = thermo.composition().Y();

    // Restore the original Ydefault if it exists, and create a new Ydefault if
    // it does not
    if (YdefaultIoExists)
    {
        fileHandler().mv
        (
            YdefaultIo.path(time_.instance())
           /"molesToMassFractions-" + YdefaultIo.name(),
            YdefaultIo.filePath()
        );
    }
    else
    {
        Ydefault.write();
    }

    // One-mole constant for conversions
    static const dimensionedScalar oneMole(dimMoles, 1);

    // Construct lists of specie molar mass, fields of specie mass, and a field
    // of total mass
    List<dimensionedScalar> W(Y.size());
    PtrList<volScalarField> m(Y.size());
    volScalarField mTotal
    (
        IOobject("mTotal", time_.timeName(), mesh_),
        mesh_,
        dimensionedScalar(dimMass, 0)
    );
    bool fromMoleFractions = false, fromMoles = false;
    forAll(Y, i)
    {
        W[i].dimensions().reset(dimMass/dimMoles);
        W[i].value() = thermo.composition().W(i);

        IOobject YIo
        (
            Y[i].name(),
            time_.timeName(),
            mesh_,
            IOobject::MUST_READ
        );
        IOobject XIo
        (
            "X_" + Y[i].name(),
            time_.timeName(),
            mesh_,
            IOobject::MUST_READ
        );
        IOobject nIo
        (
            "n_" + Y[i].name(),
            time_.timeName(),
            mesh_,
            IOobject::MUST_READ
        );

        // Check consistency of input
        if (YIo.headerOk())
        {
            FatalErrorInFunction
                << "Mass fraction field " << YIo.name()
                << " already exists on disk" << exit(FatalError);
        }
        fromMoleFractions = fromMoleFractions || XIo.headerOk();
        fromMoles = fromMoles || nIo.headerOk();
        if (fromMoleFractions && fromMoles)
        {
            FatalErrorInFunction
                << "Mole fraction fields and moles fields "
                << " both found on disk"
                << exit(FatalError);
        }

        // Sum the contributions
        if (XIo.headerOk())
        {
            m.set(i, W[i]*oneMole*volScalarField(XIo, mesh_));
        }
        if (nIo.headerOk())
        {
            m.set(i, W[i]*volScalarField(nIo, mesh_));
        }
        if (XIo.headerOk() || nIo.headerOk())
        {
            mTotal += m[i];
        }
    }

    // Steal the thermo's mass fraction fields and delete the thermo
    Y_.transfer(thermo.composition().Y());
    thermoPtr.clear();

    // Divide the specie masses by the total mass to get the mass fractions
    forAll(Y_, i)
    {
        if (m.set(i))
        {
            Y_[i].forceAssign(m[i]/mTotal);
        }
        else
        {
            Y_.set(i, nullptr);
        }
    }

    return true;
}


bool Foam::functionObjects::massFractions::write()
{
    forAll(Y_, i)
    {
        if (Y_.set(i))
        {
            Y_[i].write();
        }
    }

    return true;
}


// ************************************************************************* //
