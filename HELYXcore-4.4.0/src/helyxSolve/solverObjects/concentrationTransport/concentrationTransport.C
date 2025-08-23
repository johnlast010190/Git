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
    (c) 2011-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "concentrationTransport.H"
#include "solverObjects/solverOption/SolverOption.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "db/dictionary/dictionary.H"
#include "db/Time/Time.H"
#include "finiteVolume/fvc/fvc.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "fields/Fields/oneField/oneField.H"
#include "cfdTools/general/boundarySetup/boundarySetup.H"
#include "sets/topoSetSource/topoSetSource.H"
#include "cfdTools/general/bound/bound.H"
#include "derivedFvPatchFields/phaseChangeHumidity/phaseChangeHumidityFvPatchScalarField.H"
#include "algorithms/subCycle/subCycle.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(concentrationTransport, 0);
}
}

makeFvSolverOption(concentrationTransport);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const volScalarField& Foam::fv::concentrationTransport::phase() const
{
    if (phaseFieldPtr_)
    {
        return *phaseFieldPtr_;
    }

    FatalErrorInFunction
        << "phaseFieldPtr_ is uninitialised. Please check phaseName entry and"
        << " ensure a field of that name exists." << nl << exit(FatalError);

    return volScalarField::New("zero", mesh_, dimensionedScalar(dimless, 0));
}


tmp<volScalarField> Foam::fv::concentrationTransport::diffusivity() const
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    volScalarField& Dt =
        obr_.lookupObjectRef<volScalarField>("Dt" + fieldName_);

    tmp<volScalarField> Deff;

    if (obr_.foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const cmpTurbModel& turb =
            obr_.lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

        Dt = turb.mut()/Sct_();
        Dt.correctBoundaryConditions();

        if (Sc_.valid())
        {
            Deff = volScalarField::New("Deff", Dt + turb.mu()/Sc_());
        }
        else if (D_().dimensions() == Dt.dimensions())
        {
            Deff = volScalarField::New("Deff", Dt + D_());
        }
        else
        {
            Deff = volScalarField::New("Deff", Dt + turb.rho()*D_());
        }
    }
    else if (obr_.foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        const icoTurbModel& turb =
            obr_.lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

        Dt = turb.nut()/Sct_();
        Dt.correctBoundaryConditions();

        if (Sc_.valid())
        {
            Deff = volScalarField::New("Deff", Dt + turb.nu()/Sc_());
        }
        else
        {
            Deff = volScalarField::New("Deff", Dt + D_());
        }
    }
    else
    {
        // No turbulence model available, return warning
        WarningInFunction
            << "A valid turbulence model could not be found in the database."
            << nl << "Continue with zero diffusivity" << endl;

        Deff =
            volScalarField::New
            (
                "Deff",
                obr_,
                mesh_,
                dimensionedScalar(D_().dimensions(), 0)
            );
    }

    // Wet area ratio correction for humidity
    // do not move as diffusivity is called
    // by fvOptions thermalHumiditySource
    const volScalarField& field =
        obr_.lookupObject<volScalarField>(fieldName_);

    typedef phaseChangeHumidityFvPatchScalarField humidityBoundary;

    forAll(field.boundaryField(), patchi)
    {
        if (isA<humidityBoundary>(field.boundaryField()[patchi]))
        {
            Deff.ref().boundaryFieldRef()[patchi] *=
                dynamic_cast<const humidityBoundary&>
                (
                    field.boundaryField()[patchi]
                ).wetAreaRatio();
        }
    }

    if (noDiffusionBoundaries_.size())
    {
        // set wall diffusion to zero for conservation
        forAll(Deff.ref().boundaryField(), patchi)
        {
            if
            (
                noDiffusionBoundaries_.found
                (
                    Deff.ref().mesh().boundaryMesh()[patchi].name()
                )
            )
            {
                Deff.ref().boundaryFieldRef()[patchi] = 0.0;
            }
        }
    }

    return Deff;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::concentrationTransport::concentrationTransport
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),
    fieldName_(dict.lookup("fieldName")),
    phiName_(dict.lookupOrDefault(word("phi"), word("phi"))),
    rhoName_(dict.lookupOrDefault(word("rho"), word("rho"))),
    D_(nullptr),
    Sct_(nullptr),
    Sc_(nullptr),
    phaseName_(nullptr),
    phaseFieldPtr_(nullptr),
    solveSpeedup_(0),
    nSubCycles_(0),
    residualPhaseField_("res", dimless, -1),
    noDiffusionBoundaries_
    (
        dict.lookupOrDefault<hashedWordList>
        (
            "noDiffusionBoundaries",
            wordList()
        )
    ),
    stablePhasic_(false),
    phasicSources_(false)
{
    // TODO: zeroBoundaryDiffusion should be optimally defined by explicit
    //       diffusion boundaries. However, temporary solution could be to
    //       to introduce set of boundaries where diffusion should be zero.
    concentrationTransport::read(dict);
    Info<< endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::concentrationTransport::~concentrationTransport()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::concentrationTransport::initialise()
{
    // Read molecular diffusivity and turbulent Schmidt number
    const word subDictName(fieldName_ + type());
    if (obr_.foundObject<fluidThermo>(basicThermo::dictName))
    {
        const fluidThermo& thermo =
            obr_.lookupObject<fluidThermo>(basicThermo::dictName);

        if (basicThermo::dictName == basicThermo::matDictName)
        {
            D_.reset
            (
                new dimensionedScalar
                (
                    "D",
                    thermo.properties().optionalSubDict(subDictName).lookup("D")
                )
            );
            Sct_.reset
            (
                new dimensionedScalar
                (
                    "Sct",
                    thermo.properties().optionalSubDict(subDictName).lookup("Sct")
                )
            );
            if (thermo.properties().optionalSubDict(subDictName).found("Sc"))
            {
                Sc_.reset
                (
                    new dimensionedScalar
                    (
                        "Sc",
                        thermo.properties().optionalSubDict(subDictName).lookup("Sc")
                    )
                );
            }
        }
        else
        {
            D_.reset
            (
                new dimensionedScalar
                (
                    "D" + fieldName_,
                    thermo.properties().lookup("D" + fieldName_)
                )
            );
            Sct_.reset
            (
                new dimensionedScalar
                (
                    "Sct" + fieldName_,
                    thermo.properties().lookup("Sct" + fieldName_)
                )
            );
            if (thermo.properties().found("Sc" + fieldName_))
            {
                Sc_.reset
                (
                    new dimensionedScalar
                    (
                        "Sc" + fieldName_,
                        thermo.properties().lookup("Sc" + fieldName_)
                    )
                );
            }
        }
    }
    else if (obr_.foundObject<dictionary>("transportProperties"))
    {
        const dictionary& transport =
            obr_.lookupObject<dictionary>("transportProperties");

        D_.reset
        (
            new dimensionedScalar
            (
                "D" + fieldName_,
                transport.lookup("D" + fieldName_)
            )
        );
        Sct_.reset
        (
            new dimensionedScalar
            (
                "Sct" + fieldName_,
                transport.lookup("Sct" + fieldName_)
            )
        );
        if (transport.found("Sc" + fieldName_))
        {
            Sc_.reset
            (
                new dimensionedScalar
                (
                    "Sc" + fieldName_,
                    transport.lookup("Sc" + fieldName_)
                )
            );
        }
    }
    else
    {
        // No turbulence model available, exit with FatalError
        FatalErrorInFunction
            << "A valid turbulence model could not be found in the database."
            << abort(FatalError);
    }

    // set pointer to phase fraction field
    if (phaseName_.valid())
    {
        phaseFieldPtr_ =
            obr_.lookupObjectPtr<volScalarField>(phaseName_());
    }

    return true;
}


void Foam::fv::concentrationTransport::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    solveNames.append(fieldName_);
    // Needs phi
    requiredDependencies.insert(fieldName_, {"U"});
    optionalDependencies.insert(fieldName_, {"fvMesh"});
    correctorMembers.insert
    (
        solverObject::outerCorrectorName, solveNames
    );
}


void Foam::fv::concentrationTransport::correct
(
    const word& solveName,
    const word& regionName
)
{
    if (debug)
    {
        Info<< "    " << "Solving passive scalar transport variable: "
            << fieldName_ <<  endl;
    }

    fv::options& fvOptions = this->fvOptions();

    volScalarField& field =
        obr_.lookupObjectRef<volScalarField>(fieldName_);

    const surfaceScalarField& phi =
        obr_.lookupObject<surfaceScalarField>(phiName_);

    // modify deltaT
    Time& time = const_cast<Time&>(mesh_.time());
    scalar deltaT = time.deltaTValue();

    if (nSubCycles_ > 1)
    {
        time.setDeltaT(deltaT*solveSpeedup_);

        for
        (
            subCycle<volScalarField> scalarSubCycle(field, nSubCycles_);
            !(++scalarSubCycle).end();
        )
        {
            solveEquation(field, phi, fvOptions);
        }
    }
    else
    {
        solveEquation(field, phi, fvOptions);
    }
/*
    bound
    (
        field,
        dimensionedScalar(fieldName_, field.dimensions(), 0.0)
    );
*/
}


void Foam::fv::concentrationTransport::solveEquation
(
    volScalarField& field,
    const surfaceScalarField& phi,
    fv::options& fvOptions
)
{
    // Calculate the diffusivity
    tmp<volScalarField> D = diffusivity();

    // setup phasic treatment
    if (phaseName_.valid())
    {
        // Better use interpolate(alpha)*interpolate(D) as e.g. in
        // reactingEulerFoam MulticomponentPhaseModel ?
        if (stablePhasic_)
        {
            residualPhaseField_ = dimensionedScalar(dimless, 0);
            D.ref() *= pos(phase() - 0.99);
        }
        else
        {
            D.ref() *= phase();
        }
    }

    if (phi.dimensions() == dimVelocity * dimArea && !phaseFieldPtr_)
    {
        tmp<fvScalarMatrix> Eq
        (
            fvm::ddt(field)
          + fvm::div
            (
                phi,
                field,
                "div(phi," + Foam::word(field.name()) + ")"
            )
          - fvm::laplacian
            (
                D,
                field,
                "laplacian(Deff," + Foam::word(field.name()) + ")"
            )
         ==
            fvOptions(field)
        );

        Eq->relax();
        fvOptions.constrain(Eq.ref());

        Eq->solve();
        fvOptions.correct(field);
    }
    else //variable density
    {
        const volScalarField& rho =
            obr_.lookupObject<volScalarField>(rhoName_);

        tmp<fvScalarMatrix> Eq
        (
            (
                (!phaseFieldPtr_ || stablePhasic_)
              ? fvm::ddt(rho, field)
              : fvm::ddt(phase(), rho, field)
            )
          + fvm::div
            (
                phaseFieldPtr_ ? (fvc::interpolate(rho)*phi)() : phi,
                field,
                "div(phi," + Foam::word(field.name()) + ")"
            )
          - fvm::laplacian
            (
                D,
                field,
                "laplacian(Deff," + Foam::word(field.name()) + ")"
            )
         ==
            (
                phasicSources_
              ? phase()*fvOptions(rho, field)
              : fvOptions(rho, field)
            )
        );

        if (residualPhaseField_.value() > 0)
        {
            Eq.ref() +=
                residualPhaseField_
               *(
                    fvm::ddt(rho, field)
                  - fvc::ddt(rho, field)
                );
        }

        Eq->relax();
        fvOptions.constrain(Eq.ref());

        Eq->solve();
        fvOptions.correct(field);
    }
}


void Foam::fv::concentrationTransport::write()
{
    // transported field is creatd with auto-write so nothing to do
    Info<< endl;
}


void Foam::fv::concentrationTransport::read(const dictionary& dict)
{
    // create field if it does not already exist
    if (!obr_.foundObject<volScalarField>(fieldName_))
    {
        IOobject fieldHeader
        (
            fieldName_,
            mesh_.time().timeName(),
            obr_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        // read if present
        if (fieldHeader.typeHeaderOk<volScalarField>(true))
        {
            autoPtr<volScalarField> field
            (
                new volScalarField(fieldHeader, mesh_)
            );
            field->store(field);
        }
        else
        {
            FatalErrorInFunction
                << "Entry fieldDefinitions is no longer supported."
                << " Please setup the field and boundary conditions"
                << " the standard way."
                << nl << exit(FatalError);
        }
    }

    // Create turbulent diffusivity if it doesn't already exist
    if (!obr_.foundObject<volScalarField>("Dt" + fieldName_))
    {
        autoPtr<volScalarField> field
        (
            new volScalarField
            (
                IOobject
                (
                    "Dt" + fieldName_,
                    mesh_.time().timeName(),
                    obr_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
        field->store(field);
    }

    //check if a phase name has been specified
    if (dict.found("phaseName"))
    {
        phaseName_.reset(new word(dict.lookup("phaseName")));

        // read additional phasic settings
        stablePhasic_ = dict.lookupOrDefault<Switch>("stablePhasic", false);
        phasicSources_ = dict.lookupOrDefault<Switch>("phasicSources", false);

        if (!stablePhasic_)
        {
            residualPhaseField_ = dimensionedScalar
            (
                "residualPhaseField",
                dimless,
                dict.lookupOrDefault<scalar>("residualPhaseField", 1e-5)
            );
        }

        // reset default phiName and rhoName to VOF phasic
        phiName_ =
            dict.lookupOrDefault
            (
                word("phi"),
                IOobject::groupName("alphaPhiv", IOobject::group(phaseName_()))
            );
        rhoName_ =
            dict.lookupOrDefault
            (
                word("rho"),
                IOobject::groupName("thermo:rho", IOobject::group(phaseName_()))
            );

        // sanity check
        if (phiName_ == "phi")
        {
            WarningInFunction
                << "Please check phi entry. Phasic transport with flux"
                << " phi does not seem an appropriate choice." << endl;
        }
    }

    solveSpeedup_ = dict.lookupOrDefault<label>("solveSpeedup", 1);
    nSubCycles_ = dict.lookupOrDefault<label>("nSubCycles", solveSpeedup_);
    noDiffusionBoundaries_ =
        dict.lookupOrDefault<hashedWordList>
        (
            "noDiffusionBoundaries",
            wordList()
        );
}



// ************************************************************************* //
