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
    (c) 2019-2023 OpenFOAM Foundation

Class
    Foam::surfaceArrheniusReactionRate

Description

    A modified Arrhenius reaction rate given by:
        k = (A * T^beta * exp(-Ta/T))*a
    Where a is the surface area per unit volume, the name of which is
    specified by the user.

SourceFiles
    surfaceArrheniusReactionRateI.H

\*---------------------------------------------------------------------------*/

#ifndef surfaceArrheniusReactionRate_H
#define surfaceArrheniusReactionRate_H

#include "reaction/reactionRate/ArrheniusReactionRate/ArrheniusReactionRate.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class surfaceArrheniusReactionRate;

Ostream& operator<<(Ostream&, const surfaceArrheniusReactionRate&);

/*---------------------------------------------------------------------------*\
                Class surfaceArrheniusReactionRate Declaration
\*---------------------------------------------------------------------------*/

class surfaceArrheniusReactionRate
:
    public ArrheniusReactionRate
{
    // Private Data

        //- Name of the surface area per unit volume field
        const word AvName_;

        //- Reference to the object registry
        const objectRegistry& obr_;

        //- The surface area per unit volume field
        mutable tmp<volScalarField::Internal> tAv_;


public:

    // Constructors

        //- Construct from dictionary
        inline surfaceArrheniusReactionRate
        (
            const speciesTable& species,
            const objectRegistry& ob,
            const dictionary& dict
        );


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "surfaceArrhenius";
        }

        //- Pre-evaluation hook
        inline void preEvaluate() const;

        //- Post-evaluation hook
        inline void postEvaluate() const;

        //- Evaluate the rate
        inline scalar operator()
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;

        //- Evaluate the derivative
        inline scalar ddT
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;

        //- Write to stream
        inline void write(Ostream& os) const;


    // Ostream Operator

        inline friend Ostream& operator<<
        (
            Ostream&,
            const surfaceArrheniusReactionRate&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "surfaceArrheniusReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //