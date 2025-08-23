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

\*---------------------------------------------------------------------------*/

#include "hConst.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hConst, 0);
    addToRunTimeSelectionTable(materialModel, hConst, dictionary);
    addNamedToRunTimeSelectionTable(materialModel, hConst, dictionary, CpConst);
    addToRunTimeSelectionTable(coefficientMixture, hConst, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hConst::hConst
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
:
    materialModel(obr, dict, phaseName, specieName, name)
{
    sMod_.setSize(modelsEnumSize_);
    hConst::read();
}


Foam::hConst::hConst
(
    const UPtrList<baseModels<scalar>>& models,
    const scalarField& coeffs,
    const word& name
)
:
    materialModel
    (
        dynamic_cast<const hConst&>(models[0]).obr_,
        dictionary::null,
        dynamic_cast<const hConst&>(models[0]).phaseName_,
        word::null,
        name
    ),
    Cp_(0),
    hf_(0),
    Tref_(dynamic_cast<const hConst&>(models[0]).Tref_),
    hsref_(0),
    Sf_(0)
{
    forAll(models, i)
    {
        const hConst& comp = dynamic_cast<const hConst&>(models[i]);
        if (notEqual(Tref_, comp.Tref_))
        {
            FatalErrorInFunction
                << "Tref " << Tref_<< " for "
                << (this->name().size() ? this->name() : "others")
                << " != " << comp.Tref_ << " for "
                << (comp.name().size() ? comp.name() : "others")
                << exit(FatalError);
        }
        Cp_ += coeffs[i]*comp.Cp_;
        hf_ += coeffs[i]*comp.hf_;
        hsref_ += coeffs[i]*comp.hsref_;
        Sf_ += coeffs[i]*comp.Sf_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hConst::updateTable(const word& modelName)
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
    const List<wordList> depNames
    ({
        {CpContributionModel::typeName},
        {hContributionModel::typeName},
        {hContributionModel::typeName},
        {spContributionModel::typeName},
        {CpModel::typeName, CpMCvModel::typeName},
        {rhoModel::typeName, hsModel::typeName},
        {rhoModel::typeName, haModel::typeName}
    });

    // Dependency indices
    const List<labelList> depInds
    ({
        {CpContributionInd},
        {hContributionInd},
        {hContributionInd},
        {spContributionInd},
        {CpInd, CpMCvInd},
        {rhoInd, hsInd},
        {rhoInd, haInd},
    });

    // Create links to dependent models
    fill(modelName, modList, depNames, depInds);
}


Foam::baseModels<Foam::scalar>* Foam::hConst::castScalarModel
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

    return nullptr;
}


bool Foam::hConst::read()
{
    initialise("T");
    initialise("p");
    Cp_ = dict_->lookup<scalar>("Cp");
    hf_ = dict_->lookup<scalar>("Hf");

    // Foundation has linearisation around Tstd by default,
    // but this needs to be investigated since it is changing results
    // for non-reacting cases.
    Tref_ =
        dict_->found("Tref")
      ? dict_->lookup<scalar>("Tref") + T_->offset().value()
      : 0;

    hsref_ = dict_->lookupOrDefault<scalar>("Hsref", 0);
    Sf_ = dict_->lookupOrDefault<scalar>("Sf", 0);

    return true;
}


// ************************************************************************* //
