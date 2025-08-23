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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "clouds/derived/basicReactingMultiphaseCloud/basicReactingMultiphaseCloud.H"

#include "parcels/include/makeParcelCloudFunctionObjects.H"

// Kinematic
#include "parcels/include/makeThermoParcelForces.H" // thermo variant
#include "parcels/include/makeParcelDispersionModels.H"
#include "parcels/include/makeReactingMultiphaseParcelInjectionModels.H" // MP variant
#include "parcels/include/makeParcelPatchInteractionModels.H"
#include "parcels/include/makeReactingMultiphaseParcelStochasticCollisionModels.H" // MP variant
#include "parcels/include/makeReactingParcelSurfaceFilmModels.H" // Reacting variant

// Thermodynamic
#include "parcels/include/makeParcelHeatTransferModels.H"

// Reacting
#include "parcels/include/makeReactingMultiphaseParcelCompositionModels.H" // MP Variant
#include "parcels/include/makeReactingParcelPhaseChangeModels.H"

// Reacting multiphase
#include "parcels/include/makeReactingMultiphaseParcelDevolatilisationModels.H"
#include "parcels/include/makeReactingMultiphaseParcelSurfaceReactionModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeParcelCloudFunctionObjects(basicReactingMultiphaseCloud);

// Kinematic sub-models
makeThermoParcelForces(basicReactingMultiphaseCloud);
makeParcelDispersionModels(basicReactingMultiphaseCloud);
makeReactingMultiphaseParcelInjectionModels(basicReactingMultiphaseCloud);
makeParcelPatchInteractionModels(basicReactingMultiphaseCloud);
makeReactingMultiphaseParcelStochasticCollisionModels
(
    basicReactingMultiphaseCloud
);
makeReactingParcelSurfaceFilmModels(basicReactingMultiphaseCloud);

// Thermo sub-models
makeParcelHeatTransferModels(basicReactingMultiphaseCloud);

// Reacting sub-models
makeReactingMultiphaseParcelCompositionModels
(
    basicReactingMultiphaseCloud
);
makeReactingParcelPhaseChangeModels(basicReactingMultiphaseCloud);

// Reacting multiphase sub-models
makeReactingMultiphaseParcelDevolatilisationModels
(
    basicReactingMultiphaseCloud
);
makeReactingMultiphaseParcelSurfaceReactionModels
(
    basicReactingMultiphaseCloud
);


// ************************************************************************* //
