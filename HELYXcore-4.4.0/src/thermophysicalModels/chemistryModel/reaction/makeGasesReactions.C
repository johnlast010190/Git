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
    (c) 2012-2020 OpenFOAM Foundation
    (c) 2017-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "include/dummyThermo.H"
#include "reaction/makeReaction.H"
#include "reaction/reactionRate/ArrheniusReactionRate/ArrheniusReactionRate.H"
#include "reaction/reactionRate/LandauTellerReactionRate/LandauTellerReactionRate.H"
#include "reaction/reactionRate/thirdBodyArrheniusReactionRate/thirdBodyArrheniusReactionRate.H"
#include "reaction/reactionRate/ChemicallyActivatedReactionRate/ChemicallyActivatedReactionRate.H"
#include "reaction/reactionRate/JanevReactionRate/JanevReactionRate.H"
#include "reaction/reactionRate/powerSeries/powerSeriesReactionRate.H"
#include "reaction/reactionRate/FallOffReactionRate/FallOffReactionRate.H"
#include "reaction/reactionRate/fallOffFunctions/LindemannFallOffFunction/LindemannFallOffFunction.H"
#include "reaction/reactionRate/fallOffFunctions/SRIFallOffFunction/SRIFallOffFunction.H"
#include "reaction/reactionRate/fallOffFunctions/TroeFallOffFunction/TroeFallOffFunction.H"
#include "include/forCommonGases.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    forCommonGases(defineReaction, nullArg);

    // Irreversible/reversible/non-equilibrium-reversible reactions
    forCommonGases(makeIRNReactions, ArrheniusReactionRate);
    forCommonGases(makeIRNReactions, LandauTellerReactionRate);
    forCommonGases(makeIRNReactions, thirdBodyArrheniusReactionRate);

    // Irreversible/reversible reactions
    forCommonGases(makeIRReactions, JanevReactionRate);
    forCommonGases(makeIRReactions, powerSeriesReactionRate);

    // Pressure dependent reactions
    forCommonGases
    (
        makeIRTemplate2Reactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    );
    forCommonGases
    (
        makeIRTemplate2Reactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    );
    forCommonGases
    (
        makeIRTemplate2Reactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    );
    forCommonGases
    (
        makeIRTemplate2Reactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    );
    forCommonGases
    (
        makeIRTemplate2Reactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    );
    forCommonGases
    (
        makeIRTemplate2Reactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    );
}

// ************************************************************************* //
