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
    (c) held by original author

Class
    Foam::relaxationSchemes::relaxationSchemeSpatialImplicit

Description
    An object which results in relaxation, which is only dependent on the
    spatial location within a relaxationShape.

SourceFiles
    relaxationSchemeSpatialImplicit.C

Author
    Niels Gjøl Jacobsen, Deltares.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef relaxationSchemeSpatialImplicit_H
#define relaxationSchemeSpatialImplicit_H

#include "relaxationZone/relaxationScheme/relaxationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationSchemes
{

/*---------------------------------------------------------------------------*\
                           Class relaxationSchemeSpatialImplicit Declaration
\*---------------------------------------------------------------------------*/

class relaxationSchemeSpatialImplicit
:
    public relaxationScheme
{

protected:

    // Protected data
        scalar factor_;
        scalarField weight_;

    // Protected member functions
        volScalarField& relaxationWeightsMomentum();

       volScalarField& targetAlpha();

        volVectorField& targetVelocity();

public:

    //- Runtime type information
    TypeName("relaxationSchemeSpatialImplicit");

    // Constructors

        //- from components
        relaxationSchemeSpatialImplicit
        (
            const word& ,
            const fvMesh& mesh_,
            vectorField& U,
            scalarField& alpha
        );

    // Destructor

        virtual ~relaxationSchemeSpatialImplicit()
        {}

    // Member Functions
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationSchemes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
