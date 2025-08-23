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
    (c) 2017 OpenFOAM Foundation
    (c) 2022-2024 Engys Ltd.

Class
    Foam::Function1Types::quadraticRamp

Description
    Quadratic ramp function starting from 0 and increasing quadraticRampally
    to 1 from \c t_0 over the \c duration and remaining at 1 thereafter.

See also
    Foam::Function1Types::ramp

SourceFiles
    quadraticRamp.C

\*---------------------------------------------------------------------------*/

#ifndef quadraticRamp_H
#define quadraticRamp_H

#include "primitives/functions/Function1/ramp/ramp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class quadraticRamp Declaration
\*---------------------------------------------------------------------------*/

class quadraticRamp
:
    public ramp
{
public:

    // Runtime type information
    TypeName("quadraticRamp");


    // Constructors

        //- Construct from entry name and dictionary
        quadraticRamp
        (
            const word& entryName,
            const dictionary& dict
        );


    //- Destructor
    virtual ~quadraticRamp();


    // Member Functions

        //- Return value for time t
        scalar value(const scalar t) const;

        //- Derivative
        scalar derivative(const scalar t) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const quadraticRamp&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
