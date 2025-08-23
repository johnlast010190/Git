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
    (c) 2024-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "ePolynomial.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ePolynomial, 0);
    addToRunTimeSelectionTable(materialModel, ePolynomial, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ePolynomial::ePolynomial
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
    ePolynomial::read();
}


Foam::ePolynomial::ePolynomial
(
    const UPtrList<baseModels<scalar>>& models,
    const scalarField& coeffs,
    const word& name
)
:
    materialModel
    (
        dynamic_cast<const ePolynomial&>(models[0]).obr_,
        dictionary::null,
        dynamic_cast<const ePolynomial&>(models[0]).phaseName_,
        word::null,
        name
    ),
    hf_(0),
    Sf_(0)
{
    const ePolynomial& comp = dynamic_cast<const ePolynomial&>(models[0]);
    CvCoeffs_.reset(comp.CvCoeffs_);
    eCoeffs_.reset(comp.eCoeffs_);
    sCoeffs_.reset(comp.sCoeffs_);
    hf_ = comp.hf_;
    Sf_ = comp.Sf_;
    for (label i = 1; i < models.size(); ++i)
    {
        const ePolynomial& comp = dynamic_cast<const ePolynomial&>(models[i]);
        CvCoeffs_ += coeffs[i]*comp.CvCoeffs_;
        eCoeffs_ += coeffs[i]*comp.eCoeffs_;
        sCoeffs_ += coeffs[i]*comp.sCoeffs_;
        hf_ += coeffs[i]*comp.hf_;
        Sf_ += coeffs[i]*comp.Sf_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ePolynomial::updateTable(const word& modelName)
{
    // Model names
    const wordList modList
    ({
        CvModel::typeName,
        eaModel::typeName,
        esModel::typeName,
        CpModel::typeName,
        hsModel::typeName,
        haModel::typeName,
        sModel::typeName
    });
    dep_.setSize(modList.size());

    // Required dependencies
    const List<wordList> depNames
    ({
        {CvContributionModel::typeName},
        {eContributionModel::typeName},
        {eContributionModel::typeName},
        {CvModel::typeName, CpMCvModel::typeName},
        {esModel::typeName, rhoModel::typeName},
        {eaModel::typeName, rhoModel::typeName},
        {svContributionModel::typeName}
    });

    // Dependency indices
    const List<labelList> depInds
    ({
        {CvContributionInd},
        {eContributionInd},
        {eContributionInd},
        {CvInd, CpMCvInd},
        {esInd, rhoInd},
        {eaInd, rhoInd},
        {svContributionInd}
    });

    // Create links to dependent models
    fill(modelName, modList, depNames, depInds);
}


Foam::baseModels<Foam::scalar>* Foam::ePolynomial::castScalarModel
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


bool Foam::ePolynomial::read()
{
    initialise("T");
    initialise("p");
    hf_ = dict_->lookup<scalar>("Hf");
    Sf_ = dict_->lookup<scalar>("Sf");
    CvCoeffs_.reset(dict_->lookup<scalarField>("CpCoeffs"));
    eCoeffs_.reset(CvCoeffs_.integral());
    eCoeffs_.reset(CvCoeffs_.integralMinus1());

    // Offset s poly so that it is relative to the entropy at Tstd
    sCoeffs_[0] += Sf_ - sCoeffs_.value(Tstd);
    return true;
}


// ************************************************************************* //
