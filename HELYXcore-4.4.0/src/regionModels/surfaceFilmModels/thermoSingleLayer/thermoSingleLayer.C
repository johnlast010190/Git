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
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2022-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "surfaceFilmModels/thermoSingleLayer/thermoSingleLayer.H"
#include "finiteVolume/fvc/fvcDdt.H"
#include "finiteVolume/fvc/fvcDiv.H"
#include "finiteVolume/fvc/fvcFlux.H"
#include "finiteVolume/fvm/fvmDdt.H"
#include "finiteVolume/fvm/fvmDiv.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "fields/fvPatchFields/derived/mappedField/mappedFieldFvPatchField.H"
#include "meshes/polyMesh/polyDistributionMap/distributionMap.H"
#include "global/constants/constants.H"
#include "fields/fvPatchFields/basic/mixed/mixedFvPatchFields.H"
#include "basicThermo/basicThermo.H"
#include "mixtures/basicSpecieMixture/basicSpecieMixture.H"
#include "fluidThermo/fluidThermo.H"
#include "surfaceFilmModels/submodels/thermo/heatTransferModel/heatTransferModel/heatTransferModel.H"
#include "surfaceFilmModels/submodels/thermo/phaseChangeModel/phaseChangeModel/phaseChangeModel.H"
#include "surfaceFilmModels/submodels/thermo/filmRadiationModel/filmRadiationModel/filmRadiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermoSingleLayer, 0);
addToRunTimeSelectionTable(surfaceFilmRegionModel, thermoSingleLayer, mesh);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool thermoSingleLayer::read()
{
    // No additional properties to read
    return kinematicSingleLayer::read();
}


void thermoSingleLayer::resetPrimaryRegionSourceTerms()
{
    DebugInFunction << endl;

    kinematicSingleLayer::resetPrimaryRegionSourceTerms();

    hSpPrimary_.forceAssign(dimensionedScalar(hSp_.dimensions(), 0));
}


void thermoSingleLayer::correctHforMappedT()
{
    volScalarField& T = thermo_->T();

    T.correctBoundaryConditions();
    volScalarField::Boundary& heBf = thermo_->he().boundaryFieldRef();

    forAll(heBf, patchi)
    {
        const fvPatchField<scalar>& Tp = T.boundaryField()[patchi];
        if (isA<mappedFieldFvPatchField<scalar>>(Tp))
        {
            heBf[patchi].forceAssign(thermo().he(Tp, patchi));
        }
    }
}


void thermoSingleLayer::transferPrimaryRegionThermoFields()
{
    DebugInFunction << endl;

    kinematicSingleLayer::transferPrimaryRegionThermoFields();

    // Update primary region fields on local region via direct mapped (coupled)
    // boundary conditions
    TPrimary_.correctBoundaryConditions();
    forAll(YPrimary_, i)
    {
        YPrimary_[i].correctBoundaryConditions();
    }
}


void thermoSingleLayer::transferPrimaryRegionSourceFields()
{
    DebugInFunction << endl;

    kinematicSingleLayer::transferPrimaryRegionSourceFields();

    volScalarField::Boundary& hSpPrimaryBf = hSpPrimary_.boundaryFieldRef();

    // Convert accumulated source terms into per unit area per unit time
    const scalar deltaT = time_.deltaTValue();
    forAll(hSpPrimaryBf, patchi)
    {
        scalarField rpriMagSfdeltaT
        (
            (1/deltaT)/primaryMesh().magSf().boundaryField()[patchi]
        );

        hSpPrimaryBf[patchi] *= rpriMagSfdeltaT;
    }

    // Retrieve the source fields from the primary region
    toRegion(hSp_, hSpPrimaryBf);
    hSp_.field() /= VbyA();
}


void thermoSingleLayer::correctCoverage()
{
    if (hydrophilic_)
    {
        const scalar hydrophilicDry = hydrophilicDryScale_*deltaWet_;
        const scalar hydrophilicWet = hydrophilicWetScale_*deltaWet_;

        forAll(coverage_, i)
        {
            if ((coverage_[i] < 0.5) && (delta_[i] > hydrophilicWet))
            {
                coverage_[i] = 1;
            }
            else if ((coverage_[i] > 0.5) && (delta_[i] < hydrophilicDry))
            {
                coverage_[i] = 0;
            }
        }

        coverage_.correctBoundaryConditions();
    }
    else
    {
        coverage_.forceAssign
        (
            pos(delta_ - dimensionedScalar(dimLength, deltaWet_))
        );
    }
}


void thermoSingleLayer::updateSubmodels()
{
    DebugInFunction << endl;

    // Update heat transfer coefficient sub-models
    htcs_->correct();
    htcw_->correct();

    // Update radiation
    radiation_->correct();

    // Update ejection model - mass returned is mass available for ejection
    ejection_.correct(availableMass_, cloudMassTrans_, cloudDiameterTrans_);

    phaseChange_->correct
    (
        time_.deltaTValue(),
        availableMass_,
        primaryMassTrans_,
        primaryEnergyTrans_
    );

    const volScalarField::Internal rMagSfDt((1/time().deltaT())/magSf());

    // Vapour recoil pressure
    pSp_ -= sqr(rMagSfDt*primaryMassTrans_())/(2*rhoPrimary_());

    // Update transfer model - mass returned is mass available for transfer
    transfer_.correct
    (
        availableMass_,
        primaryMassTrans_,
        primaryMomentumTrans_,
        primaryEnergyTrans_
    );

    const volScalarField::Internal rVDt
    (
        1/(time().deltaT()*regionMesh().V())
    );

    volScalarField& he = thermo_->he();

    // Update source fields
    rhoSp_ += rVDt*(cloudMassTrans_() + primaryMassTrans_());
    USp_ += rVDt*(cloudMassTrans_()*U_() + primaryMomentumTrans_());
    hSp_ += rVDt*(cloudMassTrans_()*he() + primaryEnergyTrans_());

    turbulence_->correct();
}


tmp<fvScalarMatrix> thermoSingleLayer::q(volScalarField& he) const
{
    const volScalarField::Internal& T = thermo().T();
    const volScalarField::Internal& Cpv = thermo().Cpv();

    return
    (
        // Heat-transfer to the primary region
      - fvm::Sp((htcs_->h()/VbyA())/Cpv, he)
      + (htcs_->h()/VbyA())*(he()/Cpv + coverage_()*(TPrimary_() - T))

        // Heat-transfer to the wall
      - fvm::Sp((htcw_->h()/VbyA())/Cpv, he)
      + (htcw_->h()/VbyA())*(he()/Cpv + coverage_()*(Tw() - T))
    );
}


void thermoSingleLayer::solveEnergy()
{
    DebugInFunction << endl;

    correctHforMappedT();

    volScalarField& he = thermo_->he();

    fvScalarMatrix heEqn
    (
        fvm::ddt(alpha_, rho(), he) + fvm::div(phi_, he)
      - fvm::Sp(continuityErr_, he)
     ==
      - hSp_
      + q(he)
      + radiation_->Shs()/VbyA()
    );

    heEqn.relax();

    heEqn.solve();

    thermo_->correct();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermoSingleLayer::thermoSingleLayer
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g,
    const word& regionType,
    const bool readFields
)
:
    kinematicSingleLayer(modelType, mesh, g, regionType, false),
    phaseName_(coeffs_.lookupOrDefault("phase", word::null)),
    primaryThermo_
    (
        mesh.lookupObject<fluidThermo>
        (
            IOobject::groupName(basicThermo::dictName, phaseName_)
        )
    ),
    primaryEnergyTrans_
    (
        IOobject("primaryEnergyTrans", time().timeName(), regionMesh()),
        regionMesh(),
        dimensionedScalar(dimEnergy, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    deltaWet_(coeffs_.lookup<scalar>("deltaWet")),
    hydrophilic_(coeffs_.lookup<bool>("hydrophilic")),
    hydrophilicDryScale_(0.0),
    hydrophilicWetScale_(0.0),
    hSp_
    (
        IOobject("hSp", time().timeName(), regionMesh()),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
    ),
    hSpPrimary_
    (
        IOobject(hSp_.name(), time().timeName(), primaryMesh()),
        primaryMesh(),
        dimensionedScalar(hSp_.dimensions(), 0)
    ),
    TPrimary_
    (
        IOobject("T", time().timeName(), regionMesh()),
        regionMesh(),
        dimensionedScalar(dimTemperature, 0),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),
    YPrimary_(),
    htcs_
    (
        heatTransferModel::New(*this, coeffs().subDict("upperSurfaceModels"))
    ),
    htcw_
    (
        heatTransferModel::New(*this, coeffs().subDict("lowerSurfaceModels"))
    ),
    phaseChange_(phaseChangeModel::New(*this, coeffs())),
    radiation_(radiationModel::New(*this, coeffs())),
    Tmin_(-VGREAT),
    Tmax_(VGREAT)
{
    if (coeffs().readIfPresent("Tmin", Tmin_))
    {
        Info<< "    limiting minimum temperature to " << Tmin_ << endl;
    }

    if (coeffs().readIfPresent("Tmax", Tmax_))
    {
        Info<< "    limiting maximum temperature to " << Tmax_ << endl;
    }

    if (isA<basicSpecieMixture>(primaryThermo_))
    {
        const basicSpecieMixture& primarySpecieThermo =
            refCast<const basicSpecieMixture>(primaryThermo_);

        YPrimary_.setSize(primarySpecieThermo.species().size());

        forAll(primarySpecieThermo.species(), i)
        {
            YPrimary_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        primarySpecieThermo.species()[i],
                        time().timeName(),
                        regionMesh()
                    ),
                    regionMesh(),
                    dimensionedScalar(dimless, 0),
                    this->mappedFieldAndInternalPatchTypes<scalar>()
                )
            );
        }
    }

    if (hydrophilic_)
    {
        hydrophilicDryScale_ = coeffs_.lookup<scalar>("hydrophilicDryScale");
        hydrophilicWetScale_ = coeffs_.lookup<scalar>("hydrophilicWetScale");
    }

    if (readFields)
    {
        thermoSingleLayer::transferPrimaryRegionThermoFields();

        thermoSingleLayer::correctCoverage();

        surfaceScalarField phi
        (
            IOobject
            (
                "phi",
                time().timeName(),
                regionMesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE,
                false
            ),
            fvc::flux(alpha_*rho()*U_)
        );

        phi_.forceAssign(phi);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

thermoSingleLayer::~thermoSingleLayer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void thermoSingleLayer::addSources
(
    const label patchi,
    const label facei,
    const scalar massSource,
    const vector& momentumSource,
    const scalar pressureSource,
    const scalar energySource
)
{
    kinematicSingleLayer::addSources
    (
        patchi,
        facei,
        massSource,
        momentumSource,
        pressureSource,
        energySource
    );

    DebugInFunction << "    energy   = " << energySource << endl;

    hSpPrimary_.boundaryFieldRef()[patchi][facei] -= energySource;
}


void thermoSingleLayer::preEvolveRegion()
{
    DebugInFunction << endl;

    kinematicSingleLayer::preEvolveRegion();

    // Update phase change
    primaryEnergyTrans_.forceAssign(dimensionedScalar(dimEnergy, 0));
}


void thermoSingleLayer::evolveRegion()
{
    DebugInFunction << endl;

    // Update film coverage indicator
    correctCoverage();

    correctHforMappedT();

    // Predict "delta" from continuity
    predictDelta();

    // Update sub-models to provide updated source contributions
    updateSubmodels();

    // Predict "delta" from continuity with updated source
    predictDelta();

    // Capillary pressure
    const volScalarField pc(this->pc());

    // Temporary solution to make consistent number of pimple correctors in
    // helyxSolve
    pimple_.setCorr(0);

    while (pimple_.loop())
    {
        // External pressure
        const volScalarField pe(this->pe());

        // Solve for momentum for U
        const fvVectorMatrix UEqn(solveMomentum(pc, pe));

        // Solve energy for h also updates thermo
        solveEnergy();

        // Film thickness correction loop
        while (pimple_.correct())
        {
            solveAlpha(UEqn, pc, pe);
        }
    }

    // Reset source terms for next time integration
    resetPrimaryRegionSourceTerms();
}


tmp<volScalarField::Internal> thermoSingleLayer::Ts() const
{
    return thermo().T();
}


tmp<volScalarField::Internal> thermoSingleLayer::Tw() const
{
    tmp<volScalarField::Internal> tTw
    (
        volScalarField::Internal::New
        (
            "Tw",
            regionMesh(),
            dimensionedScalar(dimTemperature, 0)
        )
    );

    volScalarField::Internal& Tw = tTw.ref();

    const volScalarField& T = thermo().T();

    // Push boundary film temperature into wall temperature internal field
    for (label i=0; i<intCoupledPatchIDs_.size(); i++)
    {
        label patchi = intCoupledPatchIDs_[i];
        const polyPatch& pp = regionMesh().boundaryMesh()[patchi];
        UIndirectList<scalar>(Tw, pp.faceCells()) = T.boundaryField()[patchi];
    }

    return tTw;
}


void thermoSingleLayer::info()
{
    kinematicSingleLayer::info();

    const scalarField& Tinternal = thermo().T();

    Info<< indent << "min/mean/max(T)    = "
        << gMin(Tinternal) << ", "
        << gAverage(Tinternal) << ", "
        << gMax(Tinternal) << nl;

    phaseChange_->info(Info);
}


tmp<volScalarField::Internal> thermoSingleLayer::Srho(const label i) const
{
    const basicSpecieMixture& primarySpecieThermo =
        refCast<const basicSpecieMixture>(primaryThermo_);

    const word evaporationSpecies =
        thermo_->properties().lookup<word>("evaporationSpecies");
    const label vapId = primarySpecieThermo.species()[evaporationSpecies];

    tmp<volScalarField::Internal> tSrho
    (
        volScalarField::Internal::New
        (
            IOobject::modelName("Srho(" + Foam::name(i) + ")", typeName),
            primaryMesh(),
            dimensionedScalar(dimMass/dimVolume/dimTime, 0)
        )
    );

    if (vapId == i)
    {
        scalarField& Srho = tSrho.ref();
        const scalarField& V = primaryMesh().V();
        const scalar dt = time().deltaTValue();

        forAll(intCoupledPatchIDs_, i)
        {
            const label filmPatchi = intCoupledPatchIDs_[i];

            scalarField patchMass =
                primaryMassTrans_.boundaryField()[filmPatchi];

            toPrimary(filmPatchi, patchMass);

            const label primaryPatchi = primaryPatchIDs()[i];
            const labelUList& cells =
                primaryMesh().boundaryMesh()[primaryPatchi].faceCells();

            forAll(patchMass, j)
            {
                Srho[cells[j]] += patchMass[j]/(V[cells[j]]*dt);
            }
        }
    }

    return tSrho;
}


tmp<volScalarField::Internal> thermoSingleLayer::Sh() const
{
    tmp<volScalarField::Internal> tSh
    (
        volScalarField::Internal::New
        (
            IOobject::modelName("Sh", typeName),
            primaryMesh(),
            dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
        )
    );

    scalarField& Sh = tSh.ref();
    const scalarField& V = primaryMesh().V();
    const scalar dt = time_.deltaTValue();

    forAll(intCoupledPatchIDs_, i)
    {
        const label filmPatchi = intCoupledPatchIDs_[i];

        scalarField patchEnergy =
            primaryEnergyTrans_.boundaryField()[filmPatchi];

        toPrimary(filmPatchi, patchEnergy);

        const UList<label>& cells =
            primaryMesh().boundaryMesh()[primaryPatchIDs()[i]].faceCells();

        forAll(patchEnergy, j)
        {
            Sh[cells[j]] += patchEnergy[j]/(V[cells[j]]*dt);
        }
    }

    return tSh;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // end namespace Foam
} // end namespace regionModels
} // end namespace surfaceFilmModels

// ************************************************************************* //
