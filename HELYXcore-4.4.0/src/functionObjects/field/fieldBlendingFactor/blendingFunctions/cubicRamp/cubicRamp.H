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
    (c) 2011-2024 Engys Ltd.

Class
    Foam::cubicRamp

Description
    cubicRamp blending functions

SourceFiles
    cubicRamp.C

\*---------------------------------------------------------------------------*/

#ifndef cubicRamp_H
#define cubicRamp_H

#include "fieldBlendingFactor/blendingFunctions/blendingFunction/blendingFunction.H"
#include "primitives/functions/Polynomial/Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingFunctions
{

/*---------------------------------------------------------------------------*\
                          Class cubicRamp Declaration
\*---------------------------------------------------------------------------*/

class cubicRamp
:
    public blendingFunction
{
private:

    // Private Data

        //- Minimum value for blending
        scalar minValue_;

        //- Maximum value for blending
        scalar maxValue_;

        //- Blending coeff at minimum value
        scalar minCoeff_;

        //- Blending coeff at maximum value
        scalar maxCoeff_;

        //- Blending polynomial
        Polynomial<4> cubicF_;


public:

    //- Runtime type information
    TypeName("cubic");


    // Constructors

        //- Construct from components
        cubicRamp(const dictionary& data, bool dum = true);

        //- Disallow default bitwise copy construction
        cubicRamp(const cubicRamp&) = delete;


    //- Destructor
    virtual ~cubicRamp(){}


    // Member Functions

        //- Update blending function
        virtual void update(scalar);


    // Member Operators

        //- Convert value into blending factor
        virtual scalar operator()(scalar) const;

        //- Disallow default bitwise assignment
        void operator=(const cubicRamp&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blendingFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
