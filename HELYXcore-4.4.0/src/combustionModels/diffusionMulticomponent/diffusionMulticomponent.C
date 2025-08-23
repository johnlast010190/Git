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
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2022-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "diffusionMulticomponent/diffusionMulticomponent.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "mixtures/multicomponentMixture/multicomponentMixture.H"
#include "cfdTools/general/include/fvCFD.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{
    defineTypeNameAndDebug(diffusionMulticomponent, 0);
    addToRunTimeSelectionTable
    (
        combustionModel,
        diffusionMulticomponent,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

Foam::PtrList<Foam::reaction>
Foam::combustionModels::diffusionMulticomponent::loadReactions
(
    const speciesTable& species,
    const dictionary& dict
) const
{
    PtrList<reaction> reactions;
    if (dict.found("reactions"))
    {
        const dictionary& reactionsDict = dict.subDict("reactions");
        forAllConstIter(dictionary, reactionsDict, iter)
        {
            reactions.append(new reaction(species, iter().dict()));
        }
    }

    return reactions;
}


void Foam::combustionModels::diffusionMulticomponent::init()
{
    // Load default values
    this->coeffs().readIfPresent("Ci", Ci_);
    this->coeffs().readIfPresent("YoxStream", YoxStream_);
    this->coeffs().readIfPresent("YfStream", YfStream_);
    this->coeffs().readIfPresent("sigma", sigma_);
    this->coeffs().readIfPresent("ftCorr", ftCorr_);
    this->coeffs().readIfPresent("alpha", alpha_);
    this->coeffs().readIfPresent("laminarIgn", laminarIgn_);

    const speciesTable& species = this->thermo().composition().species();

    scalarList specieStoichCoeffs(species.size());
    const label nReactions = reactions_.size();

    for (label k=0; k < nReactions; k++)
    {
        RijPtr_.set
        (
            k,
            new volScalarField
            (
                IOobject
                (
                    "Rijk" + Foam::name(k),
                    this->mesh_.time().timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                this->mesh_,
                dimensionedScalar(dimMass/dimTime/dimVolume, 0),
                zeroGradientFvPatchScalarField::typeName
            )
        );

        RijPtr_[k].storePrevIter();

        const List<specieCoeffs>& lhs = reactions_[k].lhs();
        const List<specieCoeffs>& rhs = reactions_[k].rhs();

        const label fuelIndex = species[fuelNames_[k]];
        const label oxidantIndex = species[oxidantNames_[k]];

        const scalar Wu = mixture_.W(fuelIndex);
        const scalar Wox = mixture_.W(oxidantIndex);

        forAll(lhs, i)
        {
            const label specieI = lhs[i].index;
            specieStoichCoeffs[specieI] = -lhs[i].stoichCoeff;
            qFuel_[k] +=
                mixture_.hf(specieI)*lhs[i].stoichCoeff/Wu;
        }

        forAll(rhs, i)
        {
            const label specieI = rhs[i].index;
            specieStoichCoeffs[specieI] = rhs[i].stoichCoeff;
            qFuel_[k] -=
                mixture_.hf(specieI)*rhs[i].stoichCoeff/Wu;
        }

        Info<< "Fuel heat of combustion : " << qFuel_[k] << endl;

        s_[k] =
            (Wox*mag(specieStoichCoeffs[oxidantIndex]))
          / (Wu*mag(specieStoichCoeffs[fuelIndex]));

        Info<< "stoichiometric oxygen-fuel ratio : " << s_[k] << endl;

        stoicRatio_[k] = s_[k]*YfStream_[k]/YoxStream_[k];

        Info<< "stoichiometric air-fuel ratio : " << stoicRatio_[k] << endl;

        const scalar fStoich = 1.0/(1.0 + stoicRatio_[k]);

        Info<< "stoichiometric mixture fraction : " << fStoich << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::diffusionMulticomponent::diffusionMulticomponent
(
    const word& modelType,
    const fluidMulticomponentThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    combustionModel(modelType, thermo, turb, combustionProperties),
    outerCorrect_(this->coeffs().lookupOrDefault("outerCorrect", false)),
    mixture_(dynamic_cast<const basicSpecieMixture&>(this->thermo())),
    reactions_(loadReactions(mixture_.species(), this->subDict("reactions"))),
    RijPtr_(reactions_.size()),
    Ci_(reactions_.size(), 1.0),
    fuelNames_(this->coeffs().lookup("fuels")),
    oxidantNames_(this->coeffs().lookup("oxidants")),
    qFuel_(reactions_.size()),
    stoicRatio_(reactions_.size()),
    s_(reactions_.size()),
    YoxStream_(reactions_.size(), 0.23),
    YfStream_(reactions_.size(), 1.0),
    sigma_(reactions_.size(), 0.02),
    oxidantRes_(this->coeffs().lookup("oxidantRes")),
    ftCorr_(reactions_.size(), 0.0),
    alpha_(1),
    laminarIgn_(false),
    timeIndex_(-1),
    chemistryPtr_(basicChemistryModel::New(thermo))
{
    init();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModels::diffusionMulticomponent::~diffusionMulticomponent()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::combustionModels::diffusionMulticomponent::correct()
{
    if (!outerCorrect_ && timeIndex_ == this->mesh().time().timeIndex())
    {
        return;
    }

    const speciesTable& species = this->thermo().composition().species();

    const label nReactions = reactions_.size();

    PtrList<volScalarField> RijlPtr(nReactions);

    for (label reactioni=0; reactioni < nReactions; reactioni++)
    {
        RijlPtr.set
        (
            reactioni,
            new volScalarField
            (
                IOobject
                (
                    "Rijl" + word(reactioni),
                    this->mesh_.time().timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                this->mesh_,
                dimensionedScalar(dimMass/dimTime/dimVolume, 0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
        volScalarField& Rijl = RijlPtr[reactioni];

        // Obtain individual reaction rates for each reaction
        const label fuelIndex = species[fuelNames_[reactioni]];
        const PtrList<volScalarField::Internal> RR
        (
            chemistryPtr_->reactionRR(reactioni)
        );
        if (laminarIgn_)
        {
            Rijl.ref() = -RR[fuelIndex];
        }

        // Look for the fuelStoic
        const List<specieCoeffs>& rhs = reactions_[reactioni].rhs();
        const List<specieCoeffs>& lhs = reactions_[reactioni].lhs();

        // Set to zero RR's
        forAll(lhs, l)
        {
            const label lIndex = lhs[l].index;
            const_cast<volScalarField::Internal&>(chemistryPtr_->RR()[lIndex]) =
                dimensionedScalar(dimMass/dimTime/dimVolume, 0);
        }
        forAll(rhs, l)
        {
            const label rIndex = rhs[l].index;
            const_cast<volScalarField::Internal&>(chemistryPtr_->RR()[rIndex]) =
                dimensionedScalar(dimMass/dimTime/dimVolume, 0);
        }
    }

    for (label reactioni=0; reactioni < nReactions; reactioni++)
    {
        const label fuelIndex = species[fuelNames_[reactioni]];
        const label oxidantIndex = species[oxidantNames_[reactioni]];
        const volScalarField& Yfuel =
            this->thermo().composition().Y(fuelIndex);
        const volScalarField& Yox =
            this->thermo().composition().Y(oxidantIndex);
        const volScalarField ft
        (
            "ft" + Foam::name(reactioni),
            (
                s_[reactioni]*Yfuel - (Yox - YoxStream_[reactioni])
            )
           /(
                s_[reactioni]*YfStream_[reactioni] + YoxStream_[reactioni]
            )
        );
        const scalar sigma = sigma_[reactioni];
        const scalar fStoich = 1.0/(1.0 + stoicRatio_[reactioni]) + ftCorr_[reactioni];
        const volScalarField OAvailScaled
        (
            "OAvailScaled",
            Yox/max(oxidantRes_[reactioni], 1e-3)
        );
        const volScalarField preExp
        (
            "preExp" + Foam::name(reactioni),
            1.0  + sqr(OAvailScaled)
        );
        const volScalarField filter
        (
            (1.0/(sigma*sqrt(2.0*constant::mathematical::pi)))
          * exp(-sqr(ft - fStoich)/(2*sqr(sigma)))
        );
        const volScalarField topHatFilter(pos0(filter - 1e-3));
        const volScalarField prob("prob" + Foam::name(reactioni), preExp*filter);
        const volScalarField RijkDiff
        (
           "RijkDiff",
            Ci_[reactioni]*this->turbulence().muEff()*prob*
            (
                mag(fvc::grad(Yfuel) & fvc::grad(Yox))
            )
           *pos0(Yox)*pos0(Yfuel)
        );
        volScalarField& Rijk = RijPtr_[reactioni];
        if (laminarIgn_)
        {
            Rijk =
                min(RijkDiff, topHatFilter*RijlPtr[reactioni]*pos0(Yox)*pos0(Yfuel));
        }
        else
        {
            Rijk = RijkDiff;
        }
        Rijk.relax(alpha_);
        if (this->mesh_.time().outputTime() && debug)
        {
            Rijk.write();
            ft.write();
        }

        // Look for the fuelStoic
        const List<specieCoeffs>& rhs = reactions_[reactioni].rhs();
        const List<specieCoeffs>& lhs = reactions_[reactioni].lhs();
        scalar fuelStoic = 1.0;
        forAll(lhs, l)
        {
            const label lIndex = lhs[l].index;
            if (lIndex == fuelIndex)
            {
                fuelStoic = lhs[l].stoichCoeff;
                break;
            }
        }
        const scalar MwFuel = mixture_.W(fuelIndex);

        // Update left hand side species
        forAll(lhs, l)
        {
            const label lIndex = lhs[l].index;
            const scalar stoichCoeff = lhs[l].stoichCoeff;
            const_cast<volScalarField::Internal&>(chemistryPtr_->RR()[lIndex])
                += -Rijk*stoichCoeff*mixture_.W(lIndex)/fuelStoic/MwFuel;
        }

        // Update right hand side species
        forAll(rhs, r)
        {
            const label rIndex = rhs[r].index;
            const scalar stoichCoeff = rhs[r].stoichCoeff;
            const_cast<volScalarField::Internal&>(chemistryPtr_->RR()[rIndex])
                += Rijk*stoichCoeff*mixture_.W(rIndex)/fuelStoic/MwFuel;
        }
    }
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::diffusionMulticomponent::R(volScalarField& Y) const
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Y, dimMass/dimTime));

    fvScalarMatrix& Su = tSu.ref();
    const label speciei =
        this->thermo().composition().species()[Y.member()];
    Su += chemistryPtr_->RR()[speciei];

    return tSu;
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::diffusionMulticomponent::Qdot() const
{
    return volScalarField::New("Qdot", chemistryPtr_->Qdot());
}


bool Foam::combustionModels::diffusionMulticomponent::read()
{
    if (combustionModel::read())
    {
        outerCorrect_ =
            this->coeffs().lookupOrDefault("outerCorrect", false);
        this->coeffs().readIfPresent("Ci", Ci_);
        this->coeffs().readIfPresent("sigma", sigma_);
        this->coeffs().readIfPresent("oxidantRes", oxidantRes_);
        this->coeffs().readIfPresent("ftCorr", ftCorr_);
        this->coeffs().readIfPresent("alpha", alpha_);
        this->coeffs().readIfPresent("laminarIgn", laminarIgn_);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
