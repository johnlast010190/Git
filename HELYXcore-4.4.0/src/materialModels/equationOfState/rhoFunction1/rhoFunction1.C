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
    (c) 2022-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "rhoFunction1.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rhoFunction1, 0);
    addToRunTimeSelectionTable(materialModel, rhoFunction1, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rhoFunction1::rhoFunction1
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
    rhoFunction1::read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rhoFunction1::updateTable(const word& modelName)
{
    // Model names
    const wordList modList({ZModel::typeName});
    dep_.setSize(modList.size());

    // Required dependencies
    const wordList depNames(modList.size(), RModel::typeName);

    // Dependency indices
    const labelList depInds(modList.size(), RInd);

    // Create up links to dependent models
    fill(modelName, modList, depNames, depInds);
}


Foam::baseModels<Foam::scalar>* Foam::rhoFunction1::castScalarModel
(
    const word& modelName
)
{
    castMaterial(modelName, rho)
    castMaterial(modelName, hContribution)
    castMaterial(modelName, CpContribution)
    castMaterial(modelName, eContribution)
    castMaterial(modelName, CvContribution)
    castMaterial(modelName, spContribution)
    castMaterial(modelName, svContribution)
    castMaterial(modelName, psi)
    castMaterial(modelName, Z)
    castMaterial(modelName, CpMCv)
    castMaterial(modelName, alphav)
    return nullptr;
}


bool Foam::rhoFunction1::read()
{
    initialise("T");
    initialise("p");
    rhoFunction1_.clear();
    rhoFunction1_ = Function1<scalar>::New("rho", *dict_);
    rhoMin_ = dict_->lookupOrDefault<scalar>("rhoMin", 0.0);
    rhoMax_ = dict_->lookupOrDefault<scalar>("rhoMax", GREAT);

    return true;
}


// ************************************************************************* //
