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
    (c) 2022 Engys Ltd.

Class
    Foam::SoundObserver

Description

    Stores observed pressure fluctuations


SourceFiles
    SoundObserver.C

\*---------------------------------------------------------------------------*/

#ifndef soundObserver_H
#define soundObserver_H

#include "db/IOstreams/Fstreams/OFstream.H"
#include "primitives/Vector/vector/vector.H"
#include "primitives/Tuple2/Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SoundObserver Declaration
\*---------------------------------------------------------------------------*/

class SoundObserver
{
    //- Name of the sound observer
    word name_;

    //- Position of the sound observer
    vector position_;

    //- Pressure fluctuation [Pa]
    scalar pPrime_;

public:

    // Constructors

        // Construct from components
        SoundObserver
        (
            word name,
            vector pos
        );

        // Default empty object
        SoundObserver();


    // Member functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return position of observer
        const vector& position() const
        {
            return position_;
        }

        //- Return fluctuation pressure
        const scalar& pPrime() const
        {
            return pPrime_;
        }

        //- Set fluctuating pressure
        void pPrime(scalar pPrime);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
