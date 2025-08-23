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
    (c) 2016-2021 Engys Ltd.
    (c) 2011-2024 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/


#include "BirdCarreauMu.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BirdCarreauMu, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        BirdCarreauMu,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BirdCarreauMu::BirdCarreauMu
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
    BirdCarreauMu::read();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::BirdCarreauMu::updateTable(const word& modelName)
{
    // Create links to dependent models
    fill(modelName, muModel::typeName, strainRateModel::typeName, sRateInd, 0);
}


Foam::baseModels<Foam::scalar>* Foam::BirdCarreauMu::castScalarModel
(
    const word& modelName
)
{
    castMaterial(modelName, mu)

    return nullptr;
}


bool Foam::BirdCarreauMu::read()
{
    mu0_ = dict_->lookup<scalar>("mu0");
    muInf_ = dict_->lookup<scalar>("muInf");
    k_ = dict_->lookup<scalar>("k");
    n_ = dict_->lookup<scalar>("n");
    a_ = dict_->lookupOrDefault<scalar>("a", 2.0);

    return true;
}


// ************************************************************************* //
