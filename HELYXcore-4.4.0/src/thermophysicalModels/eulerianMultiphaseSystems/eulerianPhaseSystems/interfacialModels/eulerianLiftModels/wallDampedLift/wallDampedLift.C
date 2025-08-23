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
    (c) 2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "wallDampedLift.H"
#include "../../../eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianLiftModels
{
    defineTypeNameAndDebug(wallDamped, 0);
    addToRunTimeSelectionTable(eulerianLiftModel, wallDamped, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianLiftModels::wallDamped::wallDamped
(
    const dictionary& dict,
    const eulerianPhasePair& pair
)
:
    eulerianLiftModel(dict, pair),
    liftModel_(eulerianLiftModel::New(eulerianLiftModel::typeName, dict.subDict("lift"), pair)),
    wallDampingModel_
    (
        eulerianWallDampingModel::New
        (
            eulerianWallDampingModel::typeName, dict.subDict("wallDamping"), pair
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianLiftModels::wallDamped::~wallDamped()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::eulerianLiftModels::wallDamped::Cl() const
{
    return wallDampingModel_->damp(liftModel_->Cl());
}


Foam::tmp<Foam::volVectorField> Foam::eulerianLiftModels::wallDamped::Fi() const
{
    return wallDampingModel_->damp(liftModel_->Fi());
}


Foam::tmp<Foam::volVectorField> Foam::eulerianLiftModels::wallDamped::F() const
{
    return wallDampingModel_->damp(liftModel_->F());
}


Foam::tmp<Foam::surfaceScalarField> Foam::eulerianLiftModels::wallDamped::Ff() const
{
    return wallDampingModel_->damp(liftModel_->Ff());
}


// ************************************************************************* //
