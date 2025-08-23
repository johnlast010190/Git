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
    (c) 2021-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "CrossPowerLawMu.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(CrossPowerLawMu, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        CrossPowerLawMu,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CrossPowerLawMu::CrossPowerLawMu
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
    dep_.setSize(1);
    CrossPowerLawMu::read();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::CrossPowerLawMu::updateTable(const word& modelName)
{
    // Create links to dependent models
    fill(modelName, muModel::typeName, strainRateModel::typeName, sRateInd, 0);
}


Foam::baseModels<Foam::scalar>* Foam::CrossPowerLawMu::castScalarModel
(
    const word& modelName
)
{
    castMaterial(modelName, mu)

    return nullptr;
}


bool Foam::CrossPowerLawMu::read()
{
    mu0_ = dict_->lookup<scalar>("mu0");
    muInf_ = dict_->lookup<scalar>("muInf");
    m_ = dict_->lookup<scalar>("m");
    n_ = dict_->lookup<scalar>("n");

    return true;
}


// ************************************************************************* //
