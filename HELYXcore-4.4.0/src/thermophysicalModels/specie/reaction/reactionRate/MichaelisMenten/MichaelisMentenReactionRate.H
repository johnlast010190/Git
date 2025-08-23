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
    (c) 2018-2021 OpenFOAM Foundation

Class
    Foam::MichaelisMentenReactionRate

Description

    Michaelis-Menten reaction rate for enzymatic reactions.

    Reference:
    \verbatim
        Michaelis, L., & Menten, M. L. (1913).
        Die Kinetik der InwertinWirkung.
        Biochem., (49), 333-369.
    \endverbatim

SourceFiles
    MichaelisMentenReactionRateI.H

\*---------------------------------------------------------------------------*/

#ifndef MichaelisMentenReactionRate_H
#define MichaelisMentenReactionRate_H

#include "reaction/Reactions/Reaction/Reaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class MichaelisMentenReactionRate;

Ostream& operator<<(Ostream&, const MichaelisMentenReactionRate&);


/*---------------------------------------------------------------------------*\
               Class MichaelisMentenReactionRate Declaration
\*---------------------------------------------------------------------------*/

class MichaelisMentenReactionRate
{
    // Private Data

        //- List of species present in reaction system
        const speciesTable& species_;

        //- The maximum reaction rate at saturating substrate concentration
        scalar Vmax_;

        //- The Michaelis constant
        //  the substrate concentration at which the reaction rate is half Vmax_
        scalar Km_;

        //- The substrate specie index
        label s_;


public:

    // Constructors

        //- Construct from dictionary
        inline MichaelisMentenReactionRate
        (
            const speciesTable& species,
            const dictionary& dict
        );


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "MichaelisMenten";
        }

        //- Pre-evaluation hook
        inline void preEvaluate() const;

        //- Post-evaluation hook
        inline void postEvaluate() const;

        //- Return the rate
        inline scalar operator()
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;

        //- The derivative of the rate w.r.t. temperature
        inline scalar ddT
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;

        //- Is the rate a function of concentration?
        inline bool hasDdc() const;

        //- The derivative of the rate w.r.t. concentration
        inline void ddc
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li,
            scalarField& ddc
        ) const;

        //- Write to stream
        inline void write(Ostream& os) const;


    // Ostream Operator

        inline friend Ostream& operator<<
        (
            Ostream&,
            const MichaelisMentenReactionRate&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "reaction/reactionRate/MichaelisMenten/MichaelisMentenReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //