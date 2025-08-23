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
    (c) 2015-2020 OpenFOAM Foundation
    (c) 2021-2024 Engys Ltd

\*---------------------------------------------------------------------------*/
#include "humidityThermalSource.H"
#include "solverObjects/solverOption/SolverOption.H"
#include "fvMatrices/fvMatrices.H"
#include "finiteVolume/fvc/fvc.H"
#include "derivedFvPatchFields/phaseChangeHumidity/phaseChangeHumidityFvPatchScalarField.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(humidityThermalSource, 0);
}
}

makeFvSolverOption(humidityThermalSource);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::humidityThermalSource::humidityThermalSource
(
    const word& sourceName,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(sourceName, obr, dict),
    TName_("T")
{
    // If fieldNames not specified, they will attempt to be deduced from thermo
    const wordList fieldNames
    (
        coeffs_.lookupOrDefault<wordList>("fieldNames", wordList())
    );

    if (fieldNames.size() > 2)
    {
        FatalErrorInFunction
            << "Only one or two fieldNames should be specified. Settings are: "
            << fieldNames << exit(FatalError);
    }
    else if (fieldNames.size() == 2)
    {
        TName_ = fieldNames[1];
        if (TName_ == "h" || TName_ == "e")
        {
            WarningInFunction
                << "Check fieldNames order. "
                << "First one should represent volume and second boundary "
                << "source. " << nl
                << "Vomume source name: " << fieldNames.first() << nl
                << "Boundary source name: " << TName_ << nl << endl;
        }
    }

    if (fieldNames.size())
    {
        heTName_ = fieldNames.first();
    }
    transportName_ =
        dict_.lookupOrDefault<word>("humidity", "Yhumidity");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::humidityThermalSource::initialise()
{
    cTransport_ =
        obr_.lookupObjectRefPtr<concentrationTransport>(transportName_);
    thermo_ = obr_.lookupObjectPtr<fluidThermo>(basicThermo::dictName);
    transport_ =
        obr_.lookupObjectPtr<incompressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

    if (!cTransport_ || (!thermo_ && !transport_))
    {
        FatalErrorInFunction
            << "Solver object requires \"" << transportName_
            << "\" and one of: \""
            << basicThermo::dictName <<  "\", "
            << basicThermo::matDictName <<  "\", "
            << turbulenceModel::propertiesName <<  "\"."
            << exit(FatalError);
        return false;
    }

    humidityName_ = cTransport_->dict().lookup<word>("fieldName");
    w_ = obr_.lookupObjectPtr<volScalarField>(humidityName_);

    // Load humidity Cp
    word CpName("Cp");
    const word subDictName(humidityName_ + concentrationTransport::typeName);
    const dictionary* CpwDict;
    if (thermo_)
    {
        if (basicThermo::dictName != basicThermo::matDictName)
        {
            CpName += humidityName_;
        }
        CpwDict = &thermo_->properties().optionalSubDict(subDictName);
    }
    else
    {
        CpName += humidityName_;
        CpwDict = obr_.lookupObjectPtr<dictionary>("transportProperties");
        CpwDict = &CpwDict->optionalSubDict(subDictName);
    }

    Cpw_.reset(new dimensionedScalar(CpName, CpwDict->lookup(CpName)));

    return true;
}


void Foam::fv::humidityThermalSource::getSourceGraph
(
    wordList& fieldNames,
    HashTable<wordList>& sourceDependencies
)
{
    if (heTName_ == word::null)
    {
        if (thermo_)
        {
            heTName_ = thermo_->heT().name();
            TName_ = thermo_->T().name();
        }
        else
        {
            FatalIOErrorInFunction(coeffs_)
                << "No material properties object found to detect energy field "
                << "name. Please specify fieldNames."
                << exit(FatalIOError);
        }
    }


    fieldNames = {heTName_};
    sourceDependencies.insert(heTName_, {humidityName_});
}


tmp<surfaceScalarField> Foam::fv::humidityThermalSource::phiVDiff()
{
    const volScalarField Deff(cTransport_->diffusivity());

    const volScalarField& Cp =
        thermo_ ? thermo_->Cp() : transport_->transport().Cp();

    volScalarField Cprel((Cpw_() - Cp)/Cp);

    // Vapour diffusion flux
    tmp<surfaceScalarField> phiVdiff =
        fvc::interpolate(Cprel*Deff*fvc::grad(*w_)) & mesh_.Sf();

    phiVdiff.ref().rename("phiVdiff");
    return phiVdiff;
}


void Foam::fv::humidityThermalSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    //to put the term in the right side of eqn the sign must be -
    eqn -= fvm::div(this->phiVDiff(), eqn.psi());

    Info<<"fvOptions - vapor diffusion flux term added to the energy equation..."<<endl;
}


void Foam::fv::humidityThermalSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    tmp<surfaceScalarField> tPhiVDiff(this->phiVDiff());
    if (eqn.psi().dimensions() == dimEnergy/dimMass)
    {
        eqn -= fvm::div(tPhiVDiff(), thermo_->he());
    }
    else if (eqn.psi().dimensions() == dimTemperature)
    {
        eqn -=
            changeVariable
            (
                fvm::div(tPhiVDiff, thermo_->he()),
                thermo_->Cpv(),
                thermo_->T()
            );
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported dimensions for field " << eqn.psi().name() << "." << nl
            << "Must be either dimTemperature or dimEnergy/dimMass."
            << exit(FatalError);
    }

    Info<< "fvOptions - vapor diffusion flux term added to the energy equation..."
        << endl;
}


void Foam::fv::humidityThermalSource::getBoundarySourceGraph
(
    HashTable<labelList>& fieldPatchIDs,
    HashTable<wordList>& boundarySourceDependencies
)
{
    DynamicList<label> patchIDs;
    forAll(w_->boundaryField(), patchi)
    {
        if
        (
            isA<phaseChangeHumidityFvPatchScalarField>
            (
                w_->boundaryField()[patchi]
            )
        )
        {
            patchIDs.append(patchi);
        }
    }

    // This has to get added to the temperature boundary
    // even if we are solving for energy
    fieldPatchIDs.insert(TName_, patchIDs);
    boundarySourceDependencies.insert(TName_, {humidityName_});
}


void Foam::fv::humidityThermalSource::addBoundarySource
(
    const word& fieldName,
    const label patchID,
    const scalarField& pf,
    scalarField& f,
    scalarField& df
)
{
    f +=
        refCast<const phaseChangeHumidityFvPatchScalarField>
        (
            mesh_.lookupObject<volScalarField>
            (
                humidityName_
            ).boundaryField()[patchID]
        ).heatFlux();
}

// ************************************************************************* //
