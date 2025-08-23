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

#include "adiabaticPerfectFluidLaw.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(adiabaticPerfectFluidLaw, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        adiabaticPerfectFluidLaw,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adiabaticPerfectFluidLaw::adiabaticPerfectFluidLaw
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
    adiabaticPerfectFluidLaw::read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::adiabaticPerfectFluidLaw::incompressible() const
{
    return p_->isConst();
}


bool Foam::adiabaticPerfectFluidLaw::isochoric() const
{
    return p_->isConst();
}


void Foam::adiabaticPerfectFluidLaw::updateTable(const word& modelName)
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


Foam::baseModels<Foam::scalar>* Foam::adiabaticPerfectFluidLaw::castScalarModel
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


bool Foam::adiabaticPerfectFluidLaw::read()
{
    initialise("p");
    initialise("T");
    p0_ = dict_->lookup<scalar>("p0") + p_->offset().value();
    rho0_ = dict_->lookup<scalar>("rho0");
    gamma_ = dict_->lookup<scalar>("gamma");
    B_ = dict_->lookup<scalar>("B");

    return true;
}


// ************************************************************************* //
