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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2010-2015 Engys Ltd.

\*---------------------------------------------------------------------------*/

#ifndef makeParcelCollisionModels_H
#define makeParcelCollisionModels_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "submodels/Kinematic/CollisionModel/NoCollision/NoCollision.H"
#include "submodels/Kinematic/CollisionModel/PairCollision/PairCollision.H"

#include "submodels/Kinematic/CollisionModel/PairCollision/PairModel/PairSpringSliderDashpot/PairSpringSliderDashpot.H"
#include "submodels/Kinematic/CollisionModel/PairCollision/PairModel/PairJKR/PairJKR.H"
#include "submodels/Kinematic/CollisionModel/PairCollision/PairModel/PairJKRSSD/PairJKRSSD.H"

#include "submodels/Kinematic/CollisionModel/PairCollision/WallModel/WallSpringSliderDashpot/WallSpringSliderDashpot.H"
#include "submodels/Kinematic/CollisionModel/PairCollision/WallModel/WallLocalSpringSliderDashpot/WallLocalSpringSliderDashpot.H"
#include "submodels/Kinematic/CollisionModel/PairCollision/WallModel/WallLocalJKR/WallLocalJKR.H"
#include "submodels/Kinematic/CollisionModel/PairCollision/WallModel/WallLocalJKRSSD/WallLocalJKRSSD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeParcelCollisionModels(CloudType)                                   \
                                                                               \
    makeCollisionModel(CloudType);                                             \
    makeCollisionModelType(NoCollision, CloudType);                            \
    makeCollisionModelType(PairCollision, CloudType);                          \
                                                                               \
    makePairModel(CloudType);                                                  \
    makePairModelType(PairSpringSliderDashpot, CloudType);                     \
    makePairModelType(PairJKR, CloudType);                                       \
    makePairModelType(PairJKRSSD, CloudType);                                   \
                                                                               \
    makeWallModel(CloudType);                                                  \
    makeWallModelType(WallSpringSliderDashpot, CloudType);                     \
    makeWallModelType(WallLocalSpringSliderDashpot, CloudType);                \
    makeWallModelType(WallLocalJKR, CloudType);                                   \
    makeWallModelType(WallLocalJKRSSD, CloudType);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
