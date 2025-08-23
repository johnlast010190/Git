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

#include "rhoConstLaw.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rhoConstLaw, 0);
    addToRunTimeSelectionTable(materialModel, rhoConstLaw, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rhoConstLaw::rhoConstLaw
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
    rhoConstLaw::read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rhoConstLaw::updateTable(const word& modelName)
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


Foam::baseModels<Foam::scalar>* Foam::rhoConstLaw::castScalarModel
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


bool Foam::rhoConstLaw::read()
{
    initialise("p");
    rho_ = dict_->lookup<scalar>("rho");

    return true;
}


// ************************************************************************* //
