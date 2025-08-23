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
    Foam::tetherPotentials::pitchForkRing

Description


SourceFiles
    pitchForkRing.C

\*---------------------------------------------------------------------------*/

#ifndef pitchForkRing_H
#define pitchForkRing_H

#include "tetherPotential/basic/tetherPotential.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace tetherPotentials
{

/*---------------------------------------------------------------------------*\
                       Class pitchForkRing Declaration
\*---------------------------------------------------------------------------*/

class pitchForkRing
:
    public tetherPotential
{
    // Private data

        dictionary pitchForkRingCoeffs_;

        scalar mu_;
        scalar alpha_;
        scalar rOrbit_;


public:

    //- Runtime type information
    TypeName("pitchForkRing");


    // Constructors

        //- Construct from components
        pitchForkRing
        (
            const word& name,
            const dictionary& tetherPotentialProperties
        );


    //- Destructor
    ~pitchForkRing()
    {}


    // Member Functions

        scalar energy(const vector r) const;

        vector force(const vector r) const;

        //- Read dictionary
        bool read(const dictionary& tetherPotentialProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace tetherPotentials
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
