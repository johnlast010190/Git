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
    (c) 2015 OpenFOAM Foundation

Class
    Foam::eulerianWallDampingModels::sine

Description

SourceFiles
    sineWallDamping.C

\*---------------------------------------------------------------------------*/

#ifndef sineWallDamping_H
#define sineWallDamping_H

#include "../interpolated/interpolatedWallDamping.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class eulerianPhasePair;

namespace eulerianWallDampingModels
{

/*---------------------------------------------------------------------------*\
                           Class sine Declaration
\*---------------------------------------------------------------------------*/

class sine
:
    public interpolated
{
    // Private data

        //- Diameter coefficient
        const dimensionedScalar Cd_;


protected:

    // Protected member functions

        //- Return the force limiter field
        virtual tmp<volScalarField> limiter() const;


public:

    //- Runtime type information
    TypeName("sine");


    // Constructors

        //- Construct from components
        sine
        (
            const dictionary& dict,
            const eulerianPhasePair& pair
        );


    //- Destructor
    virtual ~sine();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace eulerianWallDampingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
