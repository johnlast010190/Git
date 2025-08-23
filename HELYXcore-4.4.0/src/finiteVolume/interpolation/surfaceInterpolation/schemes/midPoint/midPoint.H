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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::midPoint

Group
    grpFvSurfaceInterpolationSchemes

Description
    Mid-point interpolation (weighting factors = 0.5) scheme class.

SourceFiles
    midPoint.C

\*---------------------------------------------------------------------------*/

#ifndef midPoint_H
#define midPoint_H

#include "interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class midPoint Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class midPoint
:
    public surfaceInterpolationScheme<Type>
{
public:

    //- Runtime type information
    TypeName("midPoint");


    // Constructors

        //- Construct from mesh
        midPoint(const fvMesh& mesh)
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}

        //- Construct from Istream
        midPoint
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream&
        )
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}

        //- Construct from faceFlux and Istream
        midPoint
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream&
        )
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const VolField<Type>&
        ) const
        {
            tmp<surfaceScalarField> taw
            (
                surfaceScalarField::New
                (
                    "midPointWeights",
                    this->mesh(),
                    dimensionedScalar(dimless, 0.5)
                )
            );

            surfaceScalarField::Boundary& awBf = taw.ref().boundaryFieldRef();

            forAll(awBf, patchi)
            {
                if (!awBf[patchi].coupled())
                {
                    awBf[patchi] = 1.0;
                }
            }

            return taw;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const midPoint&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
