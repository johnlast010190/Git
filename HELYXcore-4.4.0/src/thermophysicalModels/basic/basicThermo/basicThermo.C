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
    (c) 2017 OpenCFD Ltd.
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2016-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "basicThermo/basicThermo.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "fields/fvPatchFields/derived/fixedValueZone/fixedValueZoneFvPatchField.H"
#include "fields/fvPatchFields/basic/blended/blendedFvPatchField.H"
#include "derivedFvPatchFields/fixedEnergy/fixedEnergyFvPatchScalarField.H"
#include "derivedFvPatchFields/fixedEnergyZone/fixedEnergyZoneFvPatchScalarField.H"
#include "derivedFvPatchFields/gradientEnergy/gradientEnergyFvPatchScalarField.H"
#include "derivedFvPatchFields/mixedEnergy/mixedEnergyFvPatchScalarField.H"
#include "fields/fvPatchFields/derived/mixedEnergyCalculatedTemperature/mixedEnergyCalculatedTemperatureFvPatchScalarField.H"
#include "derivedFvPatchFields/blendedEnergy/blendedEnergyFvPatchScalarField.H"
#include "fields/fvPatchFields/derived/fixedJump/fixedJumpFvPatchFields.H"
#include "fields/fvPatchFields/derived/fixedJumpAMI/fixedJumpAMIFvPatchFields.H"
#include "derivedFvPatchFields/energyJump/energyJump/energyJumpFvPatchScalarField.H"
#include "derivedFvPatchFields/energyJump/energyJumpAMI/energyJumpAMIFvPatchScalarField.H"
#include "materialModels/materialTables/materialTables.H"
#include "materialModels/baseModels/baseModels.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(basicThermo, 0);
    defineRunTimeSelectionTable(basicThermo, objectRegistry);
}

Foam::word Foam::basicThermo::dictName("thermophysicalProperties");
const Foam::word Foam::basicThermo::matDictName("materialProperties");


// * * * * * * * * * * * * Static Member Functions  * * * * * * * * * * * //

Foam::word Foam::basicThermo::heBoundaryBaseType(const fvPatchScalarField& tpf)
{
    word hbt = word::null;

    //- prioritize BCs with specified patchType
    if (tpf.patchType() != word::null)
    {
        hbt = tpf.patchType();
        return hbt;
    }

    if (tpf.overridesConstraint())
    {
        hbt = tpf.patch().type();
    }
    else if (isA<fixedJumpAMIFvPatchScalarField>(tpf))
    {
        const fixedJumpAMIFvPatchScalarField& pf =
            dynamic_cast<const fixedJumpAMIFvPatchScalarField&>(tpf);

        hbt = pf.interfaceFieldType();
    }

    return hbt;
}


Foam::word Foam::basicThermo::heBoundaryType(const fvPatchScalarField& tpf)
{
    word hbt = tpf.type();
    if (isA<fixedValueZoneFvPatchField<scalar>>(tpf))
    {
        hbt = fixedEnergyZoneFvPatchScalarField::typeName;
    }
    else if (isA<fixedValueFvPatchScalarField>(tpf))
    {
        hbt = fixedEnergyFvPatchScalarField::typeName;
    }
    else if
    (
        isA<zeroGradientFvPatchScalarField>(tpf)
     || isA<fixedGradientFvPatchScalarField>(tpf)
    )
    {
        hbt = gradientEnergyFvPatchScalarField::typeName;
    }
    else if
    (
        isA<mixedFvPatchScalarField>(tpf)
     || isA<mixedEnergyCalculatedTemperatureFvPatchScalarField>(tpf)
    )
    {
        hbt = mixedEnergyFvPatchScalarField::typeName;
    }
    else if (isA<fixedJumpFvPatchScalarField>(tpf))
    {
        hbt = energyJumpFvPatchScalarField::typeName;
    }
    else if (isA<fixedJumpAMIFvPatchScalarField>(tpf))
    {
        hbt = energyJumpAMIFvPatchScalarField::typeName;
    }
    else if (isA<blendedFvPatchField<scalar>>(tpf))
    {
        hbt = blendedEnergyFvPatchScalarField::typeName;
    }
    //else if (isA<genericFvPatchField<scalar>>(tpf))
    else if (tpf.type() == "generic")
    // so post-processing tools will work with unlinked boundaries
    {
        WarningInFunction
            << "Converting unknown fvPatchField type for patch "
            << tpf.patch().name()
            << " to `" << fixedEnergyFvPatchScalarField::typeName
            << "`" << nl
            << tab << "This is only a valid conversion for post-processing"
            << " tasks." << endl;
        hbt = fixedEnergyFvPatchScalarField::typeName;
    }
    return hbt;
}


Foam::wordList Foam::basicThermo::heBoundaryBaseTypes()
{
    const volScalarField::Boundary& tbf = T().boundaryField();
    wordList hbt(tbf.size(), word::null);
    forAll(tbf, patchi)
    {
        hbt[patchi] = heBoundaryBaseType(tbf[patchi]);
    }
    return hbt;
}


Foam::wordList Foam::basicThermo::heBoundaryTypes()
{
    const volScalarField::Boundary& tbf = T().boundaryField();
    wordList hbt(tbf.size());
    forAll(tbf, patchi)
    {
        hbt[patchi] = heBoundaryType(tbf[patchi]);
    }
    return hbt;
}


Foam::volScalarField& Foam::basicThermo::lookupOrConstruct
(
    const objectRegistry& obr,
    const word& name,
    const dimensionSet& dims,
    const scalar& defaultValue
)
{
    if (!obr.objectRegistry::foundObject<volScalarField>(name))
    {
        IOobject header
        (
            name,
            obr.time().timeName(),
            obr,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        if (header.typeHeaderOk<volScalarField>(true))
	    {
            volScalarField* fPtr(new volScalarField(header, mesh(obr)));

            // Transfer ownership of this object to the objectRegistry
            fPtr->store(fPtr);
        }
        else
        {
            volScalarField* fPtr
            (
                new volScalarField
                (
                    IOobject
                    (
                        name,
                        obr.time().timeName(),
                        obr,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(obr),
                    dimensionedScalar(name, dims, defaultValue),
                    zeroGradientFvPatchField<scalar>::typeName
                )
            );

            // Transfer ownership of this object to the objectRegistry
            fPtr->store(fPtr);
        }
    }

    return obr.lookupObjectRef<volScalarField>(name);
}


Foam::volScalarField& Foam::basicThermo::lookupOrConstruct
(
    const objectRegistry& mesh,
    const char* name,
    const dimensionSet& dims,
    const scalar& defaultValue
)
{
    return lookupOrConstruct(mesh, word(name), dims, defaultValue);
}


Foam::tmp<Foam::volScalarField> Foam::basicThermo::TAbsIfFound
(
    const objectRegistry& obr,
    const word& TName
)
{
    const basicThermo* thermoPtr =
        obr.template lookupObjectPtr<basicThermo>(basicThermo::dictName);
    return
        thermoPtr
      ? thermoPtr->TAbs()
      : volScalarField::New
        (
            "TAbs",
            obr.lookupObject<volScalarField>(TName)
        );
}


Foam::dimensionedScalar Foam::basicThermo::TRefIfFound
(
    const objectRegistry& obr
)
{
    const basicThermo* thermoPtr =
        obr.template lookupObjectPtr<basicThermo>(basicThermo::dictName);
    return
        thermoPtr
      ? thermoPtr->TRef()
      : dimensionedScalar("TRef", dimTemperature, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::volScalarField&
Foam::basicThermo::implementation::lookupOrConstructPhasic
(
    const objectRegistry& obr,
    const word& name,
    bool& found,
    const dimensionSet& dims,
    const scalar& defaultValue
)
{
    const word phasicName(IOobject::groupName(name, phaseName_));

    // Try to look up the phasic one
    found = true;

    volScalarField* Tphasic =
        obr.objectRegistry::lookupObjectRefPtr<volScalarField>(phasicName);

    if (Tphasic)
    {
        return *Tphasic;
    }

    // Try to read the phasic one
    IOobject phasicHeader
    (
        phasicName,
        obr.time().timeName(),
        obr,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (phasicHeader.typeHeaderOk<volScalarField>(true))
    {
        volScalarField* fPtr
        (
            new volScalarField(phasicHeader, mesh(obr))
        );

        // Transfer ownership of this object to the objectRegistry
        fPtr->store(fPtr);

        return *fPtr;
    }

    found = false;

    if (phaseName_ != word::null)
    {
        // Try to read the global one and initialise a 'dummy' phasic field
        // from that
        IOobject header
        (
            name,
            obr.time().timeName(),
            obr,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if (header.typeHeaderOk<volScalarField>(true))
        {
            volScalarField f0(header, mesh(obr));

            // Make a copy and turn off writing when the field did not exist in
            // setup
            phasicHeader.readOpt() = IOobject::NO_READ;
            phasicHeader.writeOpt() = IOobject::NO_WRITE;
            volScalarField* fPtr(new volScalarField(phasicHeader, f0));

            // Transfer ownership of this object to the objectRegistry
            fPtr->store(fPtr);

            return *fPtr;
        }
    }

    // Finally, create a new field; don't write
    volScalarField* fPtr
    (
        new volScalarField
        (
            IOobject
            (
                phasicName,
                obr.time().timeName(),
                obr,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(obr),
            dimensionedScalar(name, dims, defaultValue),
            zeroGradientFvPatchField<scalar>::typeName
        )
    );

    // Transfer ownership of this object to the objectRegistry
    fPtr->store(fPtr);

    return *fPtr;
}


void Foam::basicThermo::implementation::lookupAndCheckout
(
    const char* name
) const
{
    if (implementation::db().foundObject<volScalarField>(name))
    {
        implementation::db().checkOut(*implementation::db()[name]);
    }
}


Foam::IOobject Foam::basicThermo::loadDictionary
(
    const objectRegistry& obr,
    const word& phaseName,
    const bool isRegistred
)
{
    IOobject header
    (
        IOobject::groupName(dictName, phaseName),
        obr.time().caseConstant(),
        obr,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    );
    IOobject matHeader
    (
        matDictName,
        obr.time().caseSystem(),
        obr,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    );
    if
    (
        header.typeHeaderOk<dictionary>(true)
     && matHeader.typeHeaderOk<dictionary>(true)
    )
    {
        FatalErrorInFunction
            << "It isn't allowed to use both \"" << matDictName << "\""
            << " and \"thermophysicalProperties\" dictionaries."
            << nl << exit(FatalError);
    }
    if (header.typeHeaderOk<dictionary>(true))
    {
        if (dictName == matDictName)
        {
            FatalErrorInFunction
                << "It isn't allowed to use both \"" << matDictName << "\""
                << " and \"thermophysicalProperties\" dictionaries."
                << nl << exit(FatalError);
        }
        return
            IOobject
            (
                IOobject::groupName(dictName, phaseName),
                obr.time().constant(),
                obr,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                isRegistred
            );
    }
    else if (matHeader.typeHeaderOk<dictionary>(true))
    {
        dictName = matDictName;
        return
            IOobject
            (
                dictName,
                obr.time().system(),
                obr,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                isRegistred
            );
    }
    FatalErrorInFunction
        << "Either \"" << matDictName << "\""
        << " or \"thermophysicalProperties\" is required."
        << nl << exit(FatalError);
    NotImplemented;
}


Foam::basicThermo::implementation::implementation
(
    const objectRegistry& obr,
    const word& phaseName
)
:
    IOdictionary(loadDictionary(obr, phaseName, true)),
    phaseName_(phaseName),
    pRef_(readRefValue("p", dimPressure, pRefFound_)),
    TRef_(readRefValue("T", dimTemperature, TRefFound_)),
    T_(lookupOrConstructPhasic(obr, "T", TFound_, dimTemperature, 300-TRef_)),
    kappa_
    (
        IOobject
        (
            implementation::phasePropertyName("kappa", phaseName),
            obr.time().timeName(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(obr),
        dimensionedScalar(dimEnergy/dimTime/dimLength/dimTemperature, Zero)
    ),
    dpdt_(lookupOrDefault<Switch>("dpdt", true))
{
    if (dictName == matDictName)
    {
        (*this).rename(IOobject::groupName(matDictName, phaseName));
    }
}


Foam::basicThermo::implementation::implementation
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName
)
:
    IOdictionary
    (
        IOobject
        (
            IOobject::groupName(dictName, phaseName),
            obr.time().constant(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    phaseName_(phaseName),
    pRef_(readRefValue("p", dimPressure, pRefFound_)),
    TRef_(readRefValue("T", dimTemperature, TRefFound_)),
    T_(lookupOrConstructPhasic(obr, "T", TFound_, dimTemperature)),
    kappa_
    (
        IOobject
        (
            implementation::phasePropertyName("kappa", phaseName),
            obr.time().timeName(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(obr),
        dimensionedScalar(dimEnergy/dimTime/dimLength/dimTemperature, Zero)
    ),
    dpdt_(lookupOrDefault<Switch>("dpdt", true))
{
    if (dictName == matDictName)
    {
        (*this).rename(IOobject::groupName(matDictName, phaseName));
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicThermo> Foam::basicThermo::New
(
    const objectRegistry& obr,
    const word& phaseName
)
{
    return New<basicThermo>(obr, phaseName);
}


Foam::basicThermo* Foam::basicThermo::lookupPtr
(
    const objectRegistry& obr,
    const word& phaseName
)
{
    basicThermo* matPtr = obr.lookupObjectRefPtr<basicThermo>(matDictName);
    if (matPtr)
    {
        return matPtr;
    }
    else
    {
        return
            obr.lookupObjectRefPtr<basicThermo>
            (
                IOobject::groupName(dictName, phaseName)
            );
    }
}


Foam::basicThermo& Foam::basicThermo::lookupOrCreate
(
    const objectRegistry& obr,
    const word& phaseName
)
{
    basicThermo* thermoPtr = lookupPtr(obr, phaseName);
    if (!thermoPtr)
    {
        if (basicThermo::debug)
        {
            Pout<< "basicThermo::lookupOrCreate"
                << "(const objectRegistry&, const word&) : "
                << "constructing thermophysical model for region "
                << obr.name() << endl;
        }
        thermoPtr = New(obr, phaseName).ptr();
        thermoPtr->db().store
        (
            dynamic_cast<basicThermo::implementation*>(thermoPtr)
        );
    }

    return *thermoPtr;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicThermo::~basicThermo()
{}


Foam::basicThermo::implementation::~implementation()
{
    lookupAndCheckout("p");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::basicThermo::implementation::readRefValue
(
    const word& fieldName,
    const dimensionSet& dims,
    bool& refValueFound
) const
{
    refValueFound = false;
    const word refName(fieldName + "Ref");
    if (found("referenceFields") && isDict("referenceFields"))
    {
        refValueFound = subDict("referenceFields").found(fieldName);
        return
            subDict("referenceFields").lookupOrDefault<dimensionedScalar>
            (
                fieldName,
                dimensionedScalar(refName, dims, 0)
            ).value();
    }
    else if (found(refName))
    {
        refValueFound = true;
        if (!isDict("thermoType"))
        {
            Warning
                << "Loading old style " << refName << " "  << nl
                << "(not supported in new material library)"
                << nl << endl;
        }
        return lookup<scalar>(refName);
    }
    return 0.0;
}


const Foam::basicThermo& Foam::basicThermo::lookupThermo
(
    const fvPatchScalarField& pf
)
{
    word phaseDictName
    (
        IOobject::groupName(dictName, pf.internalField().group())
    );

    if (pf.db().foundObject<basicThermo>(phaseDictName))
    {
        return pf.db().lookupObject<basicThermo>(phaseDictName);
    }
    else
    {
        HashTable<const basicThermo*> thermos =
            pf.db().lookupClass<basicThermo>();

        for
        (
            HashTable<const basicThermo*>::iterator iter = thermos.begin();
            iter != thermos.end();
            ++iter
        )
        {
            if
            (
                &(iter()->he().internalField()) == &(pf.internalField())
            )
            {
                return *iter();
            }
        }
    }

    return pf.db().lookupObject<basicThermo>(phaseDictName);
}


void Foam::basicThermo::validate
(
    const string& app,
    const word& a
) const
{
    if (!(he().name() == phasePropertyName(a)))
    {
        FatalErrorInFunction
            << "Supported energy type is " << phasePropertyName(a)
            << ", thermodynamics package provides " << he().name()
            << exit(FatalError);
    }
}


void Foam::basicThermo::validate
(
    const string& app,
    const word& a,
    const word& b
) const
{
    if
    (
       !(
            he().name() == phasePropertyName(a)
         || he().name() == phasePropertyName(b)
        )
    )
    {
        FatalErrorInFunction
            << "Supported energy types are " << phasePropertyName(a)
            << " and " << phasePropertyName(b)
            << ", thermodynamics package provides " << he().name()
            << exit(FatalError);
    }
}


Foam::tmp<Foam::volScalarField> Foam::basicThermo::gamma() const
{
    return volScalarField::New(phasePropertyName("gamma"), Cp()/Cv());
}


Foam::tmp<Foam::scalarField> Foam::basicThermo::gamma
(
    const scalarField& T,
    const label patchi
) const
{
    return Cp(T, patchi)/Cv(T, patchi);
}


Foam::wordList Foam::basicThermo::splitThermoName
(
    const word& thermoName,
    const int nCmpt
)
{
    wordList cmpts(nCmpt);

    string::size_type beg=0, end=0, endb=0, endc=0;
    int i = 0;

    while
    (
        (endb = thermoName.find('<', beg)) != string::npos
     || (endc = thermoName.find(',', beg)) != string::npos
    )
    {
        if (endb == string::npos)
        {
            end = endc;
        }
        else if ((endc = thermoName.find(',', beg)) != string::npos)
        {
            end = min(endb, endc);
        }
        else
        {
            end = endb;
        }

        if (beg < end)
        {
            cmpts[i] = thermoName.substr(beg, end-beg);
            cmpts[i++].replaceAll(">","");

            // If the number of number of components in the name
            // is greater than nCmpt return an empty list
            if (i == nCmpt)
            {
                return wordList();
            }
        }
        beg = end + 1;
    }

    // If the number of number of components in the name is not equal to nCmpt
    // return an empty list
    if (i + 1 != nCmpt)
    {
        return wordList();
    }

    if (beg < thermoName.size())
    {
        cmpts[i] = thermoName.substr(beg, string::npos);
        cmpts[i].replaceAll(">","");
    }

    return cmpts;
}


const Foam::volScalarField& Foam::basicThermo::implementation::T() const
{
    return T_;
}


Foam::volScalarField& Foam::basicThermo::implementation::T()
{
    return T_;
}


const Foam::volScalarField& Foam::basicThermo::implementation::kappa() const
{
    return kappa_;
}


bool Foam::basicThermo::implementation::read()
{
    if (regIOobject::read())
    {
        pRef_ = readRefValue("p", dimPressure, pRefFound_);
        TRef_ = readRefValue("T", dimTemperature, TRefFound_);
        return true;
    }
    else
    {
        return false;
    }

}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::basicThermo::implementation::operator[](const word& modelName) const
{
    if (dictName != matDictName)
    {
        FatalErrorInFunction
            << "This operator requires material properties library"
            << exit(FatalError);
    }
    materialTables& mat =
        db().subRegistry("materialModels").lookupObjectRef<materialTables>
        (
            "materialTables"
        );

    return mat(modelName)();
}


// ************************************************************************* //
