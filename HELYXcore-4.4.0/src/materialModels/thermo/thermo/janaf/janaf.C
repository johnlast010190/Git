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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2021-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "janaf.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(janaf, 0);
    addToRunTimeSelectionTable(materialModel, janaf, dictionary);
    addToRunTimeSelectionTable(coefficientMixture, janaf, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::janaf::checkInputData() const
{
    if (Tlow_ >= Thigh_)
    {
        FatalErrorInFunction
            << "Tlow(" << Tlow_ << ") >= Thigh(" << Thigh_ << ')'
            << exit(FatalError);
    }

    if (Tcommon_ <= Tlow_)
    {
        FatalErrorInFunction
            << "Tcommon(" << Tcommon_ << ") <= Tlow(" << Tlow_ << ')'
            << exit(FatalError);
    }

    if (Tcommon_ > Thigh_)
    {
        FatalErrorInFunction
            << "Tcommon(" << Tcommon_ << ") > Thigh(" << Thigh_ << ')'
            << exit(FatalError);
    }
}


Foam::janaf::janaf
(
    const UPtrList<baseModels<scalar>>& models,
    const scalarField& coeffs,
    const word& name
)
:
    materialModel
    (
        dynamic_cast<const janaf&>(models[0]).obr_,
        dictionary::null,
        dynamic_cast<const janaf&>(models[0]).phaseName_,
        word::null,
        name
    ),
    highCpCoeffs_({0, 0, 0, 0, 0, 0, 0}),
    lowCpCoeffs_({0, 0, 0, 0, 0, 0, 0}),
    R_(GREAT),
    isConverted_(true)
{
    const janaf& comp = dynamic_cast<const janaf&>(models[0]);
    Tlow_ = comp.Tlow_;
    Thigh_ = comp.Thigh_;
    Tcommon_ = comp.Tcommon_;

    forAll(models, i)
    {
        const janaf& comp = dynamic_cast<const janaf&>(models[i]);

        Tlow_ = max(Tlow_, comp.Tlow_);
        Thigh_ = min(Thigh_, comp.Thigh_);
        if (notEqual(Tcommon_, comp.Tcommon_) && janaf::debug)
        {
            FatalErrorInFunction
                << "Tcommon " << Tcommon_ << " for "
                << (this->name().size() ? this->name() : "others")
                << " != " << comp.Tcommon_ << " for "
                << (comp.name().size() ? comp.name() : "others")
                << exit(FatalError);
        }
        forAll(highCpCoeffs_, ci)
        {
            highCpCoeffs_[ci] += coeffs[i]*comp.highCpCoeffs_[ci];
        }
        forAll(highCpCoeffs_, ci)
        {
            lowCpCoeffs_[ci] += coeffs[i]*comp.lowCpCoeffs_[ci];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::janaf::janaf
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
:
    materialModel(obr, dict, phaseName, specieName, name),
    isConverted_(false)
{
    sMod_.setSize(modelsEnumSize_);
    janaf::read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::janaf::updateTable(const word& modelName)
{
    // Model names
    const wordList modList
    ({
        CpModel::typeName,
        haModel::typeName,
        hsModel::typeName,
        sModel::typeName,
        CvModel::typeName,
        esModel::typeName,
        eaModel::typeName
    });
    dep_.setSize(modList.size());

    // Required dependencies
    const List<wordList> dependencies
    ({
        {CpContributionModel::typeName, RModel::typeName},
        {hContributionModel::typeName, RModel::typeName},
        {haModel::typeName, hfModel::typeName, RModel::typeName},
        {spContributionModel::typeName, RModel::typeName},
        {CpModel::typeName, CpMCvModel::typeName},
        {rhoModel::typeName, hsModel::typeName},
        {rhoModel::typeName, haModel::typeName}
    });

    // Dependency indices
    const List<labelList> depInds
    ({
        {CpContribution, RInd},
        {hContributionInd, RInd},
        {haInd, hfInd, RInd},
        {spContributionInd, RInd},
        {CpInd, CpMCvInd},
        {rhoInd, hsInd},
        {rhoInd, haInd},
    });

    // Create links to dependent models
    fill(modelName, modList, dependencies, depInds);

    //- Needs to be set always
    sMod_.set
    (
        RInd,
        materialTables_.sTable(phaseName_, specieName_)[RModel::typeName]
    );
    R_ = sMod_[RInd][0];

    if (!isConverted_)
    {
        isConverted_ = true;
        // Convert coefficients to mass-basis
        for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
        {
            highCpCoeffs_[coefLabel] *= R_;
            lowCpCoeffs_[coefLabel] *= R_;
        }
        checkInputData();
    }
}


Foam::baseModels<Foam::scalar>* Foam::janaf::castScalarModel
(
    const word& modelName
)
{
    castMaterial(modelName, Cp)
    castMaterial(modelName, ha)
    castMaterial(modelName, hs)
    castMaterial(modelName, hf)
    castMaterial(modelName, s)
    castMaterial(modelName, dCpdT)
    castMaterial(modelName, Cv)
    castMaterial(modelName, es)
    castMaterial(modelName, ea)
    castMaterial(modelName, limit)
    return nullptr;
}


bool Foam::janaf::read()
{
    initialise("T");
    initialise("p");
    Tlow_ = dict_->lookup<scalar>("Tlow");
    Thigh_ = dict_->lookup<scalar>("Thigh");
    Tcommon_ = dict_->lookup<scalar>("Tcommon");
    highCpCoeffs_ = dict_->lookup<scalarField>("highCpCoeffs");
    lowCpCoeffs_ = dict_->lookup<scalarField>("lowCpCoeffs");
    if (sMod_(RInd) != nullptr)
    {
        sMod_[RInd].read();
        R_ = sMod_[RInd][0];
    }

    return true;
}


// ************************************************************************* //
