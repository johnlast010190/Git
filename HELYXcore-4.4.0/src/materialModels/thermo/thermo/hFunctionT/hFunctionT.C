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
    (c) 2019-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "hFunctionT.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hFunctionT, 0);
    addToRunTimeSelectionTable(materialModel, hFunctionT, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hFunctionT::hFunctionT
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
    hFunctionT::read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hFunctionT::updateTable(const word& modelName)
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
        {haModel::typeName, hfModel::typeName},
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
        {haInd, hfInd},
        {spContributionInd},
        {CpInd, CpMCvInd},
        {rhoInd, hsInd},
        {rhoInd, haInd},
    });

    // Create links to dependent models
    fill(modelName, modList, depNames, depInds);
}


Foam::baseModels<Foam::scalar>* Foam::hFunctionT::castScalarModel
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


bool Foam::hFunctionT::read()
{
    initialise("T");
    initialise("p");
    Sf_ = dict_->lookup<scalar>("Sf");
    Cp_.clear();
    Cp_ = Function1<scalar>::New("Cp", *dict_);

    if (dict_->found("h") && !dict_->found("Hf"))
    {
        h_.clear();
        h_ = Function1<scalar>::New("h", *dict_);
        hf_ = h_->value(Tstd);
    }
    else
    {
        hf_ = dict_->lookup<scalar>("Hf");
    }

    return true;
}


// ************************************************************************* //
