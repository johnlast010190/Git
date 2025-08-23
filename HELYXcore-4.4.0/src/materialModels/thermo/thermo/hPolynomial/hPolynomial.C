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

#include "hPolynomial.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hPolynomial, 0);
    addToRunTimeSelectionTable(materialModel, hPolynomial, dictionary);
    addToRunTimeSelectionTable(coefficientMixture, hPolynomial, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hPolynomial::hPolynomial
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
    hPolynomial::read();
}


Foam::hPolynomial::hPolynomial
(
    const UPtrList<baseModels<scalar>>& models,
    const scalarField& coeffs,
    const word& name
)
:
    materialModel
    (
        dynamic_cast<const hPolynomial&>(models[0]).obr_,
        dictionary::null,
        dynamic_cast<const hPolynomial&>(models[0]).phaseName_,
        word::null,
        name
    ),
    hf_(0),
    Sf_(0)
{
    const hPolynomial& comp = dynamic_cast<const hPolynomial&>(models[0]);
    CpCoeffs_.reset(comp.CpCoeffs_);
    hCoeffs_.reset(comp.hCoeffs_);
    sCoeffs_.reset(comp.sCoeffs_);
    hf_ = comp.hf_;
    Sf_ = comp.Sf_;
    for (label i = 1; i < models.size(); ++i)
    {
        const hPolynomial& comp = dynamic_cast<const hPolynomial&>(models[i]);
        CpCoeffs_ += coeffs[i]*comp.CpCoeffs_;
        hCoeffs_ += coeffs[i]*comp.hCoeffs_;
        sCoeffs_ += coeffs[i]*comp.sCoeffs_;
        hf_ += coeffs[i]*comp.hf_;
        Sf_ += coeffs[i]*comp.Sf_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hPolynomial::updateTable(const word& modelName)
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


Foam::baseModels<Foam::scalar>* Foam::hPolynomial::castScalarModel
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


bool Foam::hPolynomial::read()
{
    initialise("T");
    initialise("p");
    hf_ = dict_->lookup<scalar>("Hf");
    Sf_ = dict_->lookup<scalar>("Sf");
    CpCoeffs_.reset(dict_->lookup<scalarField>("CpCoeffs"));
    hCoeffs_.reset(CpCoeffs_.integral());
    sCoeffs_.reset(CpCoeffs_.integralMinus1());

    // Offset s poly so that it is relative to the entropy at Tstd
    sCoeffs_[0] += Sf_ - sCoeffs_.value(Tstd);
    return true;
}


// ************************************************************************* //
