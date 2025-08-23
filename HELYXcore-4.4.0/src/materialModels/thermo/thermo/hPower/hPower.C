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
    (c) 2012-2017 OpenFOAM Foundation
    (c) 2021-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "hPower.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hPower, 0);
    addToRunTimeSelectionTable(materialModel, hPower, dictionary);
    addToRunTimeSelectionTable(coefficientMixture, hPower, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hPower::hPower
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
    hPower::read();
}


Foam::hPower::hPower
(
    const UPtrList<baseModels<scalar>>& models,
    const scalarField& coeffs,
    const word& name
)
:
    materialModel
    (
        dynamic_cast<const hPower&>(models[0]).obr_,
        dictionary::null,
        dynamic_cast<const hPower&>(models[0]).phaseName_,
        word::null,
        name
    ),
    c0_(0),
    n0_(0),
    Tref_(0),
    hf_(0)
{
    forAll(models, i)
    {
        const hPower& comp = dynamic_cast<const hPower&>(models[i]);
        c0_ += coeffs[i]*comp.c0_;
        n0_ += coeffs[i]*comp.n0_;
        Tref_ += coeffs[i]*comp.Tref_;
        hf_ += coeffs[i]*comp.hf_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hPower::updateTable(const word& modelName)
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
        {hsModel::typeName, hfModel::typeName},
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
        {hsInd, hfInd},
        {hContributionInd},
        {spContributionInd},
        {CpInd, CpMCvInd},
        {rhoInd, hsInd},
        {rhoInd, haInd},
    });

    // Create links to dependent models
    fill(modelName, modList, depNames, depInds);
}


Foam::baseModels<Foam::scalar>* Foam::hPower::castScalarModel
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


bool Foam::hPower::read()
{
    initialise("T");
    initialise("p");
    c0_ = dict_->lookup<scalar>("C0");
    n0_ = dict_->lookup<scalar>("n0");
    Tref_ = dict_->lookup<scalar>("Tref") + T_->offset().value();
    hf_ = dict_->lookup<scalar>("Hf");

    return true;
}


// ************************************************************************* //
