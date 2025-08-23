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

#include "eConst.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(eConst, 0);
    addToRunTimeSelectionTable(materialModel, eConst, dictionary);
    addNamedToRunTimeSelectionTable(materialModel, eConst, dictionary, CvConst);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eConst::eConst
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
    eConst::read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::eConst::updateTable(const word& modelName)
{
    // Model names
    const wordList models
    ({
        CpModel::typeName,
        hsModel::typeName,
        haModel::typeName,
        CvModel::typeName,
        esModel::typeName,
        eaModel::typeName,
        sModel::typeName
    });

    dep_.setSize(models.size());

    // Required dependencies
    const List<wordList> dependencies
    ({
        {CvModel::typeName, CpMCvModel::typeName},
        {esModel::typeName, rhoModel::typeName},
        {eaModel::typeName, rhoModel::typeName},
        {CvContributionModel::typeName},
        {eContributionModel::typeName},
        {eContributionModel::typeName},
        {svContributionModel::typeName, CpModel::typeName}
    });

    // Dependency indices
    const List<labelList> depInds
    ({
        {CvInd, CpMCvInd},
        {esInd, rhoInd},
        {eaInd, rhoInd},
        {CvContributionInd},
        {eContributionInd},
        {eContributionInd},
        {svContributionInd, CpInd}
    });

    // Create links to dependent models
    fill(modelName, models, dependencies, depInds);
}


Foam::baseModels<Foam::scalar>* Foam::eConst::castScalarModel
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


bool Foam::eConst::read()
{
    initialise("p");
    initialise("T");
    Cv_ = dict_->lookup<scalar>("Cv");
    hf_ = dict_->lookup<scalar>("Hf");

    // Foundation has linearisation around Tstd by default,
    // but this needs to be investigated since it is changing results
    // for non-reacting cases.
    Tref_ =
        dict_->found("Tref")
      ? dict_->lookup<scalar>("Tref") + T_->offset().value()
      : 0;

    esref_ = dict_->lookupOrDefault<scalar>("Esref", 0);

    return true;
}


// ************************************************************************* //
