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
    (c) 2011-2021 OpenFOAM Foundation

Class
    Foam::LindemannFallOffFunction

Description
    Lindemann fall-off function

SourceFiles
    LindemannFallOffFunctionI.H

\*---------------------------------------------------------------------------*/

#ifndef LindemannFallOffFunction_H
#define LindemannFallOffFunction_H

#include "primitives/Scalar/scalar/scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class LindemannFallOffFunction;
Ostream& operator<<(Ostream&, const LindemannFallOffFunction&);


/*---------------------------------------------------------------------------*\
                 Class LindemannFallOffFunction Declaration
\*---------------------------------------------------------------------------*/

class LindemannFallOffFunction
{

public:

    // Constructors

        //- Construct null
        inline LindemannFallOffFunction();

        //- Construct from dictionary
        inline LindemannFallOffFunction(const dictionary& dict);


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "Lindemann";
        }

        inline scalar operator()
        (
            const scalar T,
            const scalar Pr
        ) const;

        inline scalar ddT
        (
            const scalar T,
            const scalar Pr,
            const scalar F
        ) const;

        inline scalar ddPr
        (
            const scalar T,
            const scalar Pr,
            const scalar F
        ) const;

        //- Write to stream
        inline void write(Ostream& os) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const LindemannFallOffFunction&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "reaction/reactionRate/fallOffFunctions/LindemannFallOffFunction/LindemannFallOffFunctionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
