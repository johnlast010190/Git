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
    (c) 2011 OpenFOAM Foundation

Class
    Foam::patchInteractionData

Description
    Helper class for the LocalInteraction patch interaction model

\*---------------------------------------------------------------------------*/

#ifndef patchInteractionData_H
#define patchInteractionData_H

#include "db/IOstreams/IOstreams/Istream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class patchInteractionData Declaration
\*---------------------------------------------------------------------------*/

// Forward declaration of classes
class patchInteractionData;

// Forward declaration of friend functions
Istream& operator>>
(
    Istream& is,
    patchInteractionData& pid
);


class patchInteractionData
{
    // Private data

        //- Interaction type name
        word interactionTypeName_;

        //- Patch name
        word patchName_;

        //- Elasticity coefficient
        scalar e_;

        //- Restitution coefficient
        scalar mu_;


public:

    // Constructor

        //- Construct null
        patchInteractionData();


    // Member functions

        // Access

            //- Return const access to the interaction type name
            const word& interactionTypeName() const;

            //- Return const access to the patch name
            const word& patchName() const;

            //- Return const access to the elasticity coefficient
            scalar e() const;

            //- Return const access to the restitution coefficient
            scalar mu() const;


        // I-O

            //- Istream operator
            friend Istream& operator>>
            (
                Istream& is,
                patchInteractionData& pid
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
