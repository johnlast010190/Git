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
    (c) 2017-2021 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "surfaceFilmModels/submodels/kinematic/transferModels/transferModel/transferModel.H"
#include "surfaceFilmModels/thermoSingleLayer/thermoSingleLayer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(transferModel, 0);
defineRunTimeSelectionTable(transferModel, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void transferModel::addToTransferredMass(const scalar dMass)
{
    transferredMass_ += dMass;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

transferModel::transferModel(surfaceFilmRegionModel& film)
:
    filmSubModelBase(film),
    transferredMass_(0.0)
{}


transferModel::transferModel
(
    const word& modelType,
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    filmSubModelBase(film, dict, typeName, modelType),
    transferredMass_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

transferModel::~transferModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void transferModel::correct()
{
    if (writeTime())
    {
        scalar transferredMass0 = getModelProperty<scalar>("transferredMass");
        transferredMass0 += returnReduce(transferredMass_, sumOp<scalar>());
        setModelProperty<scalar>("transferredMass", transferredMass0);
        transferredMass_ = 0.0;
    }
}


scalar transferModel::transferredMassTotal() const
{
    scalar transferredMass0 = getModelProperty<scalar>("transferredMass");
    return transferredMass0 + returnReduce(transferredMass_, sumOp<scalar>());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
