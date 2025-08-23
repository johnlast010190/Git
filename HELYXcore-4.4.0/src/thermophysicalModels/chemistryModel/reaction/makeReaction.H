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
    (c) 2011-2023 OpenFOAM Foundation

InClass
    Foam::makeReaction

Description
    Macros for instantiating reactions on given thermo packages

\*---------------------------------------------------------------------------*/

#ifndef makeReaction_H
#define makeReaction_H

#include "reaction/Reactions/Reaction/Reaction.H"
#include "reaction/Reactions/IrreversibleReaction/IrreversibleReaction.H"
#include "reaction/Reactions/ReversibleReaction/ReversibleReaction.H"
#include "reaction/Reactions/NonEquilibriumReversibleReaction/NonEquilibriumReversibleReaction.H"
#include "thermo/thermo/thermo.H"
#include "transport/sutherland/sutherlandTransport.H"
#include "thermo/janaf/janafThermo.H"
#include "equationOfState/perfectGas/perfectGas.H"
#include "transport/polynomial/polynomialTransport.H"
#include "thermo/hPolynomial/hPolynomialThermo.H"
#include "equationOfState/icoPolynomial/icoPolynomial.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineReaction(nullArg, ThermoPhysics)                                 \
                                                                               \
    typedef Reaction<ThermoPhysics> Reaction##ThermoPhysics;                   \
                                                                               \
    defineTemplateTypeNameAndDebug(Reaction##ThermoPhysics, 0);                \
                                                                               \
    defineTemplateRunTimeSelectionTable(Reaction##ThermoPhysics, dictionary)   \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        Reaction##ThermoPhysics,                                               \
        objectRegistry                                                         \
    )


#define makeReaction(ReactionType, ReactionRate, ThermoPhysics)                \
                                                                               \
    typedef Reaction<ThermoPhysics> Reaction##ThermoPhysics;                   \
    typedef ReactionType<ThermoPhysics, ReactionRate>                          \
        ReactionType##ThermoPhysics##ReactionRate;                             \
                                                                               \
    template<>                                                                 \
    const word ReactionType##ThermoPhysics##ReactionRate::typeName             \
    (                                                                          \
        ReactionType::typeName_()                                              \
      + ReactionRate::type().capitalise()                                      \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        Reaction##ThermoPhysics,                                               \
        ReactionType##ThermoPhysics##ReactionRate,                             \
        dictionary                                                             \
    )


#define makeGeneralReaction(ReactionType, ReactionRate, ThermoPhysics)         \
                                                                               \
    typedef Reaction<ThermoPhysics> Reaction##ThermoPhysics;                   \
    typedef ReactionType<ThermoPhysics, ReactionRate>                          \
        ReactionType##ThermoPhysics##ReactionRate;                             \
                                                                               \
    template<>                                                                 \
    const word ReactionType##ThermoPhysics##ReactionRate::typeName             \
    (                                                                          \
        ReactionType::typeName_()                                              \
      + ReactionRate::type().capitalise()                                      \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        Reaction##ThermoPhysics,                                               \
        ReactionType##ThermoPhysics##ReactionRate,                             \
        objectRegistry                                                         \
    )


#define makeIReactions(ReactionRate, ThermoPhysics)                            \
                                                                               \
    makeReaction(IrreversibleReaction, ReactionRate, ThermoPhysics)


#define makeIRReactions(ReactionRate, ThermoPhysics)                           \
                                                                               \
    makeIReactions(ReactionRate, ThermoPhysics);                               \
    makeReaction(ReversibleReaction, ReactionRate, ThermoPhysics)


#define makeIRNReactions(ReactionRate, ThermoPhysics)                          \
                                                                               \
    makeIRReactions(ReactionRate, ThermoPhysics);                              \
    makeReaction(NonEquilibriumReversibleReaction, ReactionRate, ThermoPhysics)


#define makeITemplate1Reactions(ReactionRate, Param, ThermoPhysics)            \
                                                                               \
    typedef ReactionRate<Param> ReactionRate##Param;                           \
                                                                               \
    makeReaction                                                               \
    (                                                                          \
        IrreversibleReaction,                                                  \
        ReactionRate##Param,                                                   \
        ThermoPhysics                                                          \
    )


#define makeIRTemplate1Reactions(ReactionRate, Param, ThermoPhysics)           \
                                                                               \
    makeITemplate1Reactions(ReactionRate, Param, ThermoPhysics);               \
                                                                               \
    makeReaction                                                               \
    (                                                                          \
        ReversibleReaction,                                                    \
        ReactionRate##Param,                                                   \
        ThermoPhysics                                                          \
    )


#define makeITemplate2Reactions(ReactionRate, Param1, Param2, ThermoPhysics)   \
                                                                               \
    typedef ReactionRate<Param1, Param2> ReactionRate##Param1##Param2;         \
                                                                               \
    makeReaction                                                               \
    (                                                                          \
        IrreversibleReaction,                                                  \
        ReactionRate##Param1##Param2,                                          \
        ThermoPhysics                                                          \
    )


#define makeIRTemplate2Reactions(ReactionRate, Param1, Param2, ThermoPhysics)  \
                                                                               \
    makeITemplate2Reactions(ReactionRate, Param1, Param2, ThermoPhysics);      \
                                                                               \
    makeReaction                                                               \
    (                                                                          \
        ReversibleReaction,                                                    \
        ReactionRate##Param1##Param2,                                          \
        ThermoPhysics                                                          \
    )

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
