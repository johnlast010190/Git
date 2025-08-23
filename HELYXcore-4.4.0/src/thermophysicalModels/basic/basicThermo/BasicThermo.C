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
    (c) 2021-2025 Engys Ltd.
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2011-2023 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "BasicThermo.H"
#include "derivedFvPatchFields/gradientEnergy/gradientEnergyFvPatchScalarField.H"
#include "derivedFvPatchFields/mixedEnergy/mixedEnergyFvPatchScalarField.H"
#include "derivedFvPatchFields/blendedEnergy/blendedEnergyFvPatchScalarField.H"

#include "thermo/energy/matAbsoluteEnthalpy/matAbsoluteEnthalpy.H"
#include "thermo/energy/matAbsoluteInternalEnergy/matAbsoluteInternalEnergy.H"
#include "thermo/energy/matSensibleEnthalpy/matSensibleEnthalpy.H"
#include "thermo/energy/matSensibleInternalEnergy/matSensibleInternalEnergy.H"
#include "transport/strainRateLaw/strainRateLaw.H"
#include "general/scalarModelFunction1/scalarModelFunction1.H"
#include "transport/vapourDiffusivity/vapourDiffusivity.H"
#include "transport/pvInvertFunc/pvInvertFunc.H"


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicThermoType, class MixtureType>
void Foam::BasicThermo<BasicThermoType, MixtureType>::heBoundaryCorrection
(
    fvPatchScalarField& pf
)
{
    if (isA<gradientEnergyFvPatchScalarField>(pf))
    {
        refCast<gradientEnergyFvPatchScalarField>(pf).gradient() =
            pf.fvPatchField::snGrad();
    }
    else if (isA<mixedEnergyFvPatchScalarField>(pf))
    {
        refCast<mixedEnergyFvPatchScalarField>(pf).refGrad() =
            pf.fvPatchField::snGrad();
    }
    else if (isA<blendedEnergyFvPatchScalarField>(pf))
    {
        blendedEnergyFvPatchScalarField& bepf =
            refCast<blendedEnergyFvPatchScalarField>(pf);
        heBoundaryCorrection(bepf.boundaryOne());
        heBoundaryCorrection(bepf.boundaryTwo());
    }
}


template<class BasicThermoType, class MixtureType>
void Foam::BasicThermo<BasicThermoType, MixtureType>::heBoundaryCorrection
(
    volScalarField& h
)
{
    volScalarField::Boundary& hBf = h.boundaryFieldRef();

    forAll(hBf, patchi)
    {
        heBoundaryCorrection(hBf[patchi]);
    }
}


template<class BasicThermoType, class MixtureType>
Foam::materialTables&
Foam::BasicThermo<BasicThermoType, MixtureType>::matLookupOrConstruct
(
    const objectRegistry& obr,
    const dictionary& dict
)
{
    if (!obr.foundObject<objectRegistry>("materialModels"))
    {
        obr.store(new materialTables(obr, dict));
    }
    else if
    (
        !obr.subRegistry("materialModels").foundObject<materialTables>
        (
            "materialTables"
        )
    )
    {
        obr.store(new materialTables(obr, dict));
    }

    return
        obr.subRegistry
        (
            "materialModels"
        ).lookupObjectRef<materialTables>("materialTables");
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicThermoType, class MixtureType>
void Foam::BasicThermo<BasicThermoType, MixtureType>::reportFlags() const
{
    Info<< nl << "Solver settings from material properties:" << nl;

    // Isochoric?
    if (materials_.isochoric(this->phaseName_))
    {
        Info<< "\"isochoric\"";
    }
    else
    {
        Info<< "\"non-isochoric\"";
    }

    // Compressible?
    if (materials_.incompressible(this->phaseName_))
    {
        Info<< " \"incompressible\"";
    }
    else
    {
        Info<< " \"compressible\"";
    }

    // Constant Cpv?
    if (materials_.isCpvConst(this->phaseName_))
    {
        Info<< " \"constant ";
    }
    else
    {
        Info<< " \"varying ";
    }
    if (energyName(*this) == "h")
    {
        Info<< "Cp\"";
    }
    else if (energyName(*this) == "e")
    {
        Info<< "Cv\"";
    }
    else
    {
        Info<< "Cpv\"";
    }

    // Isotropic?
    if (materials_.isotropic(this->phaseName_))
    {
        Info<< " \"isotropic\"";
    }
    else
    {
        Info<< " \"non-isotropic\"";
    }

    // Buoyant?
    if (isBuoyant_)
    {
        if (isDistinctBuoyancy_)
        {
            Info<< " \"distinct buoyancy\"";
        }
        else
        {
            Info<< " \"buoyant\"";
        }
    }
    Info<< nl << endl;
}


template<class BasicThermoType, class MixtureType>
void Foam::BasicThermo<BasicThermoType, MixtureType>::initBuoyancyFlags
(
    const word& buoyantName
)
{
    matScalarTable& mat =
        materials_.sTable(this->phaseName_, word::null);
    if (mat.found("buoyancy"))
    {
        isBuoyant_ = true;
        if (mat["rho"] != mat["buoyancy"])
        {
            isDistinctBuoyancy_ = true;
        }
    }
}


template<class BasicThermoType, class MixtureType>
void Foam::BasicThermo<BasicThermoType, MixtureType>::loadLiquidModels()
{
    const wordList liqModels
    ({
        pvModel::typeName,
        hlModel::typeName,
        CpgModel::typeName,
        BModel::typeName,
        mugModel::typeName,
        kappagModel::typeName,
        sigmaModel::typeName
    });
    forAll(liqModels, modi)
    {
        const word& modelName = liqModels[modi];
        if (materials_.isModelInDict(modelName + "Model", this->phaseName_))
        {
            materials_.addSpeciesAndSpeciesMixtures
            (
                modelName + "Model",
                wordList({modelName}),
                this->phaseName_,
                modelName,
                muModel::typeName + "Model",
                wordList({modelName})
            );
        }
    }
    if (materials_.isModelInDict(DModel::typeName + "Model", this->phaseName_))
    {
        materials_.addSpeciesAndSpeciesMixtures
        (
            DModel::typeName + "Model",
            wordList({DModel::typeName}),
            this->phaseName_,
            vapourDiffusivity::typeName,
            muModel::typeName + "Model",
            wordList({DModel::typeName})
        );
    }
    if
    (
        materials_.isModelInDict
        (
            pvInvertModel::typeName + "Model",
            this->phaseName_
        )
    )
    {
        materials_.addSpeciesAndSpeciesMixtures
        (
            pvInvertModel::typeName + "Model",
            wordList({pvInvertModel::typeName}),
            this->phaseName_,
            pvInvertFunc::typeName,
            muModel::typeName + "Model",
            wordList({pvInvertModel::typeName})
        );
    }
}


template<class BasicThermoType, class MixtureType>
void Foam::BasicThermo<BasicThermoType, MixtureType>::init()
{
    if (this->TFound_ && this->phaseName_ != word::null)
    {
        // Phasic T was found on disk/in memory; record that we are using
        // phasic T
        addPhasicVariable(this->T_.member());
    }

    materials_.addSpeciesAndSpeciesMixtures<specieModels>
    (
        this->phaseName_,
        materialsSpecie::typeName,
        thermodynamics::typeName
    );
    materials_.addSpeciesAndSpeciesMixtures<equationOfState>(this->phaseName_);
    materials_.addSpeciesAndSpeciesMixtures<thermodynamics>(this->phaseName_);
    materials_.addSpeciesAndSpeciesMixtures<energyConversion>
    (
        this->phaseName_,
        standardThermo::typeName,
        thermodynamics::typeName
    );

    // Add energy entry on the specie level (without dict entry)
    // (thermodynamics has to be available)
    materials_.addSpecies
    (
        energy::typeName,
        energy::listModels(),
        this->phaseName_,
        dictionary::lookup<word>("energy"),
        thermodynamics::typeName
    );

    // Add energy entry on the phase level (without dict entry)
    // (thermodynamics has to be available)
    if (this->phaseName_ != word::null)
    {
        materials_.addGroup
        (
            {energy::typeName},
            energy::listModels(),
            this->phaseName_,
            word::null,
            dictionary::lookup<word>("energy"),
            false,
            word::null,
            wordList::null()
        );
    }

    const dictionary& phaseDict =
        this->dictionary::optionalSubDict(this->phaseName_);
    if (word(phaseDict.lookup("materialType")) == "solid")
    {
        materials_.addSpeciesAndSpeciesMixtures
        (
            kappaModel::typeName + "Model",
            wordList({kappaModel::typeName, vKappaModel::typeName}),
            this->phaseName_,
            word::null,
            word::null,
            wordList({kappaModel::typeName, vKappaModel::typeName})
        );
        if (materials_.isModelInDict(c0Model::typeName + "Model", this->phaseName_))
        {
            materials_.addSpeciesAndSpeciesMixtures
            (
                c0Model::typeName + "Model",
                wordList({c0Model::typeName}),
                this->phaseName_,
                word::null,
                word::null,
                wordList({c0Model::typeName})
            );
        }
    }
    else
    {
        if (he_.db().template foundObject<volVectorField>("U"))
        {
            // Strain rate model (requires U)
            materials_.addSpeciesAndSpeciesMixtures
            (
                strainRateModel::typeName + "Model",
                wordList({strainRateModel::typeName}),
                this->phaseName_,
                strainRateLaw::typeName,
                muModel::typeName + "Model",
                wordList({strainRateModel::typeName})
            );
        }
        else
        {
            Info<< "Material library didn't find velocity field."
                << " Non-newtonian models can't be used." << endl;
        }

        materials_.addSpeciesAndSpeciesMixtures
        (
            muModel::typeName + "Model",
            wordList({muModel::typeName}),
            this->phaseName_,
            word::null,
            word::null,
            wordList({muModel::typeName})
        );
        materials_.addSpeciesAndSpeciesMixtures
        (
            kappaModel::typeName + "Model",
            wordList({kappaModel::typeName}),
            this->phaseName_,
            word::null,
            word::null,
            wordList({kappaModel::typeName})
        );
        if (materials_.isModelInDict(c0Model::typeName + "Model", this->phaseName_))
        {
            materials_.addSpeciesAndSpeciesMixtures
            (
                c0Model::typeName + "Model",
                wordList({c0Model::typeName}),
                this->phaseName_,
                word::null,
                word::null,
                wordList({c0Model::typeName})
            );
        }
        loadLiquidModels();
    }

    materials_.addSpeciesAndSpeciesMixtures<energy>(this->phaseName_);

    if (phaseDict.found(limitModel::typeName + "Model"))
    {
        materials_.addGroup
        (
            {limitModel::typeName + "Model"},
            {limitModel::typeName},
            this->phaseName_
        );
    }

    const word buoyantName("buoyancyModel");
    const bool buoyancyModelPresent
    (
        materials_.isModelInDict(buoyantName, this->phaseName_)
    );

    // Seach for distinct buoyant model only if the case is buoyant
    // and with buoyancy rho different from rhoModel
    if (buoyancyModelPresent)
    {
        materials_.addSpeciesAndSpeciesMixtures
        (
            buoyantName,
            wordList({"buoyancy"}),
            this->phaseName_,
            word::null,
            word::null,
            wordList({rhoModel::typeName})
        );
    }

    materials_.addGroup
    (
        {TModel::typeName + "Model"},
        {
            TModel::typeName,
            ThsModel::typeName,
            ThaModel::typeName,
            TesModel::typeName,
            TeaModel::typeName
        },
        this->phaseName_,
        word::null,
        energyInverseTemperature::typeName,
        false,
        word::null,
        wordList::null()
    );

    // Update model pointers and dependency lists
    materials_.linkModelsAndUpdateTables();

    // Linking models are added only on updateTables steps
    if (buoyancyModelPresent)
    {
        initBuoyancyFlags(buoyantName);
    }

    // Check for circular model dependencies
    materials_.checkDepenencies();

    if (this->phaseName_ == word::null)
    {
        // Output material information
        materials_.reportModelsInfo();

        // Report solver flags
        reportFlags();
    }

    // Initialize all the times of he
    initHe(this->p_, he_);
}


template<class BasicThermoType, class MixtureType>
void Foam::BasicThermo<BasicThermoType, MixtureType>::initHe
(
    const volScalarField& p,
    volScalarField& he
)
{
    he.forceAssign(materials_(heModel::typeName, this->phaseName_)());
    this->heBoundaryCorrection(he);

    // Note: T does not have oldTime
    if (p.nOldTimes() > 0)
    {
        materials_.updateScalarField
        (
            BasicThermoType::phasePropertyName("p"),
            p.oldTime()
        );
        initHe(p.oldTime(), he.oldTime());
    }

    // Field always reset back to original
    if (this->p_.nOldTimes() > 0 && p.nOldTimes() == 0)
    {
        materials_.updateScalarField
        (
            BasicThermoType::phasePropertyName("p"),
            this->p_
        );
    }
}


template<class BasicThermoType, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<BasicThermoType, MixtureType>::tmpBoundaryField
(
    const scalarField& T,
    const label patchi,
    const word& modelName
) const
{
    scalarField Told;
    const bool isTsame = (&T == &this->T().boundaryField()[patchi]);
    if (!isTsame)
    {
        Told = this->T().boundaryField()[patchi];
        const_cast<volScalarField&>
        (
            this->T()
        ).boundaryFieldRef()[patchi].forceAssign(T);
    }

    tmp<scalarField> tField
    (
        materials_(modelName, this->phaseName_).boundaryField()[patchi]
    );

    // Reset T back to the original
    if (!isTsame)
    {
        const_cast<volScalarField&>
        (
            this->T()
        ).boundaryFieldRef()[patchi].forceAssign(Told);
    }

    return tField;
}


template<class BasicThermoType, class MixtureType>
Foam::word Foam::BasicThermo<BasicThermoType, MixtureType>::energyName
(
    const dictionary& dict
) const
{
    const word eName = dict.lookup<word>("energy");

    if (eName == matSensibleInternalEnergy::typeName)
    {
        return matSensibleInternalEnergy::name();
    }
    else if (eName == matSensibleEnthalpy::typeName)
    {
        return matSensibleEnthalpy::name();
    }
    else if (eName == matAbsoluteInternalEnergy::typeName)
    {
        return matAbsoluteInternalEnergy::name();
    }
    else if (eName == matAbsoluteEnthalpy::typeName)
    {
        return matAbsoluteEnthalpy::name();
    }
    else
    {
        FatalErrorInFunction
            << "Unknown energy type: " << eName << nl
            << exit(FatalError);
    }

    return word::null;
}


template<class BasicThermoType, class MixtureType>
bool Foam::BasicThermo<BasicThermoType, MixtureType>::isEnthalpy() const
{
    const word heType(energyName(this->properties()));
    if
    (
        heType == matAbsoluteEnthalpy::name()
     || heType == matSensibleEnthalpy::name()
    )
    {
        return true;
    }
    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicThermoType, class MixtureType>
Foam::BasicThermo<BasicThermoType, MixtureType>::BasicThermo
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName
)
:
    BasicThermoType(obr, phaseName),
    MixtureType(this->properties(), obr, phaseName),
    materials_(matLookupOrConstruct(obr, dict)),
    he_
    (
        (
            addPhasicVariable(energyName(this->properties())),
            volScalarField
            (
                IOobject
                (
                    phasePropertyName
                    (
                        energyName(this->properties()),
                        phaseName
                    ),
                    obr.time().timeName(),
                    obr,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    true
                ),
                this->mesh(obr),
                dimEnergy/dimMass,
                this->heBoundaryTypes(),
                this->heBoundaryBaseTypes()
            )
        )
    ),
    Cp_
    (
        IOobject
        (
            BasicThermoType::phasePropertyName("thermo-Cp", phaseName),
            obr.time().timeName(),
            obr
        ),
        this->mesh(obr),
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero)
    ),
    Cv_
    (
        IOobject
        (
            BasicThermoType::phasePropertyName("thermo-Cv", phaseName),
            obr.time().timeName(),
            obr
        ),
        this->mesh(obr),
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero)
    ),
    isBuoyant_(false),
    isDistinctBuoyancy_(false),
    calculateTFromhe_(true)
{
    init();
    BasicThermoType::init();
}


template<class BasicThermoType, class MixtureType>
Foam::BasicThermo<BasicThermoType, MixtureType>::BasicThermo
(
    const objectRegistry& obr,
    const word& phaseName
)
:
    BasicThermoType(obr, phaseName),
    MixtureType(this->phaseDict(), obr, phaseName),
    materials_(matLookupOrConstruct(obr, *this)),
    he_
    (
        (
            addPhasicVariable(energyName(this->properties())),
            volScalarField
            (
                IOobject
                (
                    phasePropertyName
                    (
                        energyName(this->properties()),
                        phaseName
                    ),
                    obr.time().timeName(),
                    obr,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    true
                ),
                this->mesh(obr),
                dimEnergy/dimMass,
                this->heBoundaryTypes(),
                this->heBoundaryBaseTypes()
            )
        )
    ),
    Cp_
    (
        IOobject
        (
            BasicThermoType::phasePropertyName("thermo-Cp", phaseName),
            obr.time().timeName(),
            obr
        ),
        this->mesh(obr),
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero)
    ),
    Cv_
    (
        IOobject
        (
            BasicThermoType::phasePropertyName("thermo-Cv", phaseName),
            obr.time().timeName(),
            obr
        ),
        this->mesh(obr),
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero)
    ),
    isBuoyant_(false),
    isDistinctBuoyancy_(false),
    calculateTFromhe_(true)
{
    init();
    BasicThermoType::init();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicThermoType, class MixtureType>
Foam::BasicThermo<BasicThermoType, MixtureType>::~BasicThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermoType, class MixtureType>
bool Foam::BasicThermo<BasicThermoType, MixtureType>::incompressible() const
{
    return materials_.incompressible(this->phaseName_);
}


template<class BasicThermoType, class MixtureType>
bool Foam::BasicThermo<BasicThermoType, MixtureType>::isochoric() const
{
    return materials_.isochoric(this->phaseName_);
}


template<class BasicThermoType, class MixtureType>
bool Foam::BasicThermo<BasicThermoType, MixtureType>::isCpvConst() const
{
    return materials_.isCpvConst(this->phaseName_);
}


template<class BasicThermoType, class MixtureType>
bool Foam::BasicThermo<BasicThermoType, MixtureType>::isotropic() const
{
    return materials_.isotropic(this->phaseName_);
}


template<class BasicThermoType, class MixtureType>
bool Foam::BasicThermo<BasicThermoType, MixtureType>::isConst() const
{
    wordList models =
    {
        muModel::typeName,
        rhoModel::typeName,
        psiModel::typeName,
        kappaModel::typeName,
        CpModel::typeName
    };

    bool isConst(true);
    forAll(models, i)
    {
        if (!materials_(models[i], this->phaseName_).isConst())
        {
            isConst = false;
        }
    }

    return isConst;
}


template<class BasicThermoType, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::BasicThermo<BasicThermoType, MixtureType>::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    // Update field pointers in the models:
    materials_.updateScalarField(BasicThermoType::phasePropertyName("p"), p);
    materials_.updateScalarField(BasicThermoType::phasePropertyName("T"), T);

    tmp<volScalarField> the(materials_(heModel::typeName, this->phaseName_)());

    // Reset model field pointers back
    materials_.updateScalarField
    (
        BasicThermoType::phasePropertyName("p"),
        this->p_
    );
    materials_.updateScalarField
    (
        BasicThermoType::phasePropertyName("T"),
        this->T_
    );

    return the;
}


template<class BasicThermoType, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<BasicThermoType, MixtureType>::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    const scalarField pOld(this->p_, cells);
    const scalarField TOld(this->T_, cells);
    UIndirectList<scalar>(this->p_, cells) = p;
    UIndirectList<scalar>(this->T_, cells) = T;

    tmp<scalarField> tCellsHE(new scalarField(cells.size()));
    scalarField& cellsHE = tCellsHE.ref();
    const baseModels<scalar>& heMatModel =
        materials_(heModel::typeName, this->phaseName_);
    forAll(cells, i)
    {
        const label celli = cells[i];
        cellsHE[i] = heMatModel[celli];
    }

    // Reset T/p back to the original
    UIndirectList<scalar>(this->T_, cells) = TOld;
    UIndirectList<scalar>(this->p_, cells) = pOld;

    return tCellsHE;
}


template<class BasicThermoType, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<BasicThermoType, MixtureType>::he
(
    const scalarField& T,
    const label patchi
) const
{
    return tmpBoundaryField(T, patchi, heModel::typeName);
}


template<class BasicThermoType, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::BasicThermo<BasicThermoType, MixtureType>::ha() const
{
    return materials_(haModel::typeName, this->phaseName_)();
}


template<class BasicThermoType, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::BasicThermo<BasicThermoType, MixtureType>::hs() const
{
    return materials_(hsModel::typeName, this->phaseName_)();
}


template<class BasicThermoType, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::BasicThermo<BasicThermoType, MixtureType>::hs
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    // Update field pointers in the models:
    materials_.updateScalarField(BasicThermoType::phasePropertyName("p"), p);
    materials_.updateScalarField(BasicThermoType::phasePropertyName("T"), T);

    tmp<volScalarField> ths(materials_(hsModel::typeName, this->phaseName_)());

    // Reset model field pointers back
    materials_.updateScalarField
    (
        BasicThermoType::phasePropertyName("p"),
        this->p_
    );
    materials_.updateScalarField
    (
        BasicThermoType::phasePropertyName("T"),
        this->T_
    );

    return ths;
}


template<class BasicThermoType, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<BasicThermoType, MixtureType>::hs
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    const scalarField pOld(this->p_, cells);
    const scalarField TOld(this->T_, cells);
    UIndirectList<scalar>(this->p_, cells) = p;
    UIndirectList<scalar>(this->T_, cells) = T;

    tmp<scalarField> tCellsHs(new scalarField(cells.size()));
    scalarField& cellsHs = tCellsHs.ref();
    const baseModels<scalar>& hsMatModel =
        materials_(hsModel::typeName, this->phaseName_);
    forAll(cells, i)
    {
        const label celli = cells[i];
        cellsHs[i] = hsMatModel[celli];
    }

    // Reset T/p back to the original
    UIndirectList<scalar>(this->T_, cells) = TOld;
    UIndirectList<scalar>(this->p_, cells) = pOld;

    return tCellsHs;
}


template<class BasicThermoType, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<BasicThermoType, MixtureType>::hs
(
    const scalarField& T,
    const label patchi
) const
{
    return tmpBoundaryField(T, patchi, hsModel::typeName);
}


template<class BasicThermoType, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::BasicThermo<BasicThermoType, MixtureType>::ha
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    // Update field pointers in the models:
    materials_.updateScalarField(BasicThermoType::phasePropertyName("p"), p);
    materials_.updateScalarField(BasicThermoType::phasePropertyName("T"), T);

    tmp<volScalarField> tha(materials_(haModel::typeName, this->phaseName_)());

    // Reset model field pointers back
    materials_.updateScalarField
    (
        BasicThermoType::phasePropertyName("p"),
        this->p_
    );
    materials_.updateScalarField
    (
        BasicThermoType::phasePropertyName("T"),
        this->T_
    );

    return tha;
}


template<class BasicThermoType, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<BasicThermoType, MixtureType>::ha
(
    const scalarField& T,
    const label patchi
) const
{
    return tmpBoundaryField(T, patchi, haModel::typeName);
}


template<class BasicThermoType, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<BasicThermoType, MixtureType>::The
(
    const scalarField& he,
    const scalarField& T0,
    const label patchi
) const
{
    scalarField Told, heOld;
    const bool isTSame = (&T0 == &this->T_.boundaryField()[patchi]);
    const bool isHeSame = (&he == &he_.boundaryField()[patchi]);
    if (!isTSame)
    {
        Told = this->T_.boundaryField()[patchi];
        this->T_.boundaryFieldRef()[patchi].forceAssign(T0);
    }
    if (!isHeSame)
    {
        heOld = he_.boundaryField()[patchi];
        const_cast<volScalarField&>(he_).boundaryFieldRef()[patchi].forceAssign
        (
            he
        );
    }

    tmp<scalarField> tField
    (
        materials_(TModel::typeName, this->phaseName_).boundaryField()[patchi]
    );

    // Reset T/he back to the original
    if (!isTSame)
    {
        this->T_.boundaryFieldRef()[patchi].forceAssign(Told);
    }
    if (!isHeSame)
    {
        const_cast<volScalarField&>(he_).boundaryFieldRef()[patchi].forceAssign
        (
            heOld
        );
    }

    return tField;
}


template<class BasicThermoType, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::BasicThermo<BasicThermoType, MixtureType>::hf() const
{
    return volScalarField::New
    (
        IOobject::groupName("hf", this->phaseName_),
        materials_(hfModel::typeName, this->phaseName_)()
       *materials_(WModel::typeName, this->phaseName_)()
    );
}


template<class BasicThermoType, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<BasicThermoType, MixtureType>::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    return tmpBoundaryField(T, patchi, CpModel::typeName);
}


template<class BasicThermoType, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<BasicThermoType, MixtureType>::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    return tmpBoundaryField(T, patchi, CvModel::typeName);
}


template<class BasicThermoType, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<BasicThermoType, MixtureType>::Cpv
(
    const scalarField& T,
    const label patchi
) const
{
    return tmpBoundaryField(T, patchi, CpvModel::typeName);
}


template<class BasicThermoType, class MixtureType>
const Foam::volScalarField&
Foam::BasicThermo<BasicThermoType, MixtureType>::Cpv() const
{
    if (isEnthalpy())
    {
        return Cp_;
    }
    else
    {
        return Cv_;
    }
}


template<class BasicThermoType, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::BasicThermo<BasicThermoType, MixtureType>::kappaEff
(
    const volScalarField& alphat
) const
{
    return volScalarField::New("kappaEff", this->kappa_ + Cp_*alphat);
}


template<class BasicThermoType, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<BasicThermoType, MixtureType>::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        this->kappa_.boundaryField()[patchi]
      + Cp_.boundaryField()[patchi]*alphat;
}


template<class BasicThermoType, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::BasicThermo<BasicThermoType, MixtureType>::alphaEff
(
    const volScalarField& alphat
) const
{
    if (isEnthalpy())
    {
        return volScalarField::New("alphaEff", this->kappa_/Cp_ + alphat);
    }
    else
    {
        return
            volScalarField::New("alphaEff", (this->kappa_ + Cp_*alphat)/Cv_);
    }
}


template<class BasicThermoType, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<BasicThermoType, MixtureType>::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    if (isEnthalpy())
    {
        return
            this->kappa_.boundaryField()[patchi]/Cp_.boundaryField()[patchi]
          + alphat;
    }
    else
    {
        return
            (
                this->kappa_.boundaryField()[patchi]
              + Cp_.boundaryField()[patchi]*alphat
            )/Cv_.boundaryField()[patchi];
    }
}


template<class BasicThermoType, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::BasicThermo<BasicThermoType, MixtureType>::buoyantRho() const
{
    return materials_("buoyancy", this->phaseName_)();
}


template<class BasicThermoType, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::BasicThermo<BasicThermoType, MixtureType>::W() const
{
    return materials_(WModel::typeName, this->phaseName_)();
}


template<class BasicThermoType, class MixtureType>
bool Foam::BasicThermo<BasicThermoType, MixtureType>::read()
{
    if (BasicThermoType::read())
    {
        materials_.updateSubDictsPtrs();
        return materials_.read();
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
