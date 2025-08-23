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
    (c) 2020 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

Class
    Foam::reaction

Description
    Reaction base-class holding the specie names and coefficients

SourceFiles
    reactionI.H
    reaction.C

\*---------------------------------------------------------------------------*/

#ifndef reaction_H
#define reaction_H

#include "speciesTable/speciesTable.H"
#include "reaction/specieCoeffs/specieCoeffs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class reaction;

inline Ostream& operator<<(Ostream&, const reaction&);

/*---------------------------------------------------------------------------*\
                          Class reaction Declaration
\*---------------------------------------------------------------------------*/

class reaction
{

public:

    // Static data

        //- Number of un-named reactions
        static label nUnNamedReactions;


private:

    // Private Data

        //- Name of reaction
        const word name_;

        //- List of specie names present in reaction system
        const speciesTable& species_;

        //- List of gas species present in reaction system
        speciesTable pyrolisisGases_;

        //- Specie info for the left-hand-side of the reaction
        List<specieCoeffs> lhs_;

        //- Specie info for the right-hand-side of the reaction
        List<specieCoeffs> rhs_;

        //- Gas specie index for the left-hand-side of the reaction
        List<specieCoeffs> glhs_;

        //- Gas specie index for the right-hand-side of the reaction
        List<specieCoeffs> grhs_;

        //- Does the reation has gaseous species
        bool gaseousSpcecies_;


    // Private Member Functions

        //- Return new reaction ID for un-named reactions
        label getNewReactionID();


public:

    // Constructors

        //- Construct from components
        reaction
        (
            const speciesTable& species,
            const List<specieCoeffs>& lhs,
            const List<specieCoeffs>& rhs
        );

        //- Construct as copy given new speciesTable
        reaction(const reaction&, const speciesTable& species);

        //- Construct from dictionary
        reaction
        (
            const speciesTable& species,
            const dictionary& dict
        );


    //- Destructor
    ~reaction()
    {}


    // Member Functions

        // Access

            //- Return the name of the reaction
            inline const word& name() const;

            //- Retun if the reaction has gaseous species
            inline bool hasGaseousSpecies() const;

            //- Return the components of the left hand side
            inline const List<specieCoeffs>& lhs() const;

            //- Return the components of the right hand side
            inline const List<specieCoeffs>& rhs() const;

            //- Access to the gas components of the left hand side
            inline const List<specieCoeffs>& grhs() const;

            //- Access to the gas components of the right hand side
            inline const List<specieCoeffs>& glhs() const;

            //- Return the specie list
            inline const speciesTable& species() const;

            //- Access to gas specie list
            inline const speciesTable& gasSpecies() const;


        // Solid Chemistry Functions

            //- Return string representation of the left of the reaction
            void solidReactionStrLeft(OStringStream&) const;

            //- Return string representation of the right of the reaction
            void solidReactionStrRight(OStringStream&) const;

            //- Return string representation of reaction
            string solidReactionStr(OStringStream&) const;


        //- Write
        void write(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const reaction&) = delete;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const reaction&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "reactionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //