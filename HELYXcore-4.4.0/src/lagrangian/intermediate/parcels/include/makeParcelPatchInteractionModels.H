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

\*---------------------------------------------------------------------------*/

#ifndef makeParcelPatchInteractionModels_H
#define makeParcelPatchInteractionModels_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "submodels/Kinematic/PatchInteractionModel/LocalInteraction/LocalInteraction.H"
#include "submodels/Kinematic/PatchInteractionModel/NoInteraction/NoInteraction.H"
#include "submodels/Kinematic/PatchInteractionModel/Rebound/Rebound.H"
#include "submodels/Kinematic/PatchInteractionModel/StandardWallInteraction/StandardWallInteraction.H"
#include "submodels/Kinematic/PatchInteractionModel/MultiInteraction/MultiInteraction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeParcelPatchInteractionModels(CloudType)                            \
                                                                               \
    makePatchInteractionModel(CloudType);                                      \
                                                                               \
    makePatchInteractionModelType(LocalInteraction, CloudType);                \
    makePatchInteractionModelType(NoInteraction, CloudType);                   \
    makePatchInteractionModelType(Rebound, CloudType);                         \
    makePatchInteractionModelType(StandardWallInteraction, CloudType);         \
    makePatchInteractionModelType(MultiInteraction, CloudType);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
