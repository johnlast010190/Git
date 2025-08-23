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
    Foam::heaviside

Description
    Heaviside blending functions

SourceFiles
    heaviside.C

\*---------------------------------------------------------------------------*/

#ifndef heaviside_H
#define heaviside_H

#include "fieldBlendingFactor/blendingFunctions/blendingFunction/blendingFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingFunctions
{

/*---------------------------------------------------------------------------*\
                          Class heaviside Declaration
\*---------------------------------------------------------------------------*/

class heaviside
:
    public blendingFunction
{
private:

    // Private Data

        //- blending coefficient for positive values
        scalar posCoeff_;

        //- blending coefficient for negative values
        scalar negCoeff_;

        //- offset applied to incoming data, default=0
        scalar intercept_;


public:

    //- Runtime type information
    TypeName("heaviside");


    // Constructors

        //- Construct from components
        heaviside(const dictionary& data, bool dum = true);

        //- Disallow default bitwise copy construction
        heaviside(const heaviside&) = delete;


    //- Destructor
    virtual ~heaviside(){}


    // Member Functions

        //- Update blending function
        virtual void update(scalar){};


    // Member Operators

        //- convert value into blending factor
        virtual scalar operator()(scalar) const;

        //- Disallow default bitwise assignment
        void operator=(const heaviside&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blendingFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
