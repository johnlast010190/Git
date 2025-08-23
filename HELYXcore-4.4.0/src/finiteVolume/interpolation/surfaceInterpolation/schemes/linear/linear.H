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
    Foam::linear

Group
    grpFvSurfaceInterpolationSchemes

Description
    Central-differencing interpolation scheme class

SourceFiles
    linear.C

\*---------------------------------------------------------------------------*/

#ifndef linear_H
#define linear_H

#include "interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class linear Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class linear
:
    public surfaceInterpolationScheme<Type>
{
public:

    //- Runtime type information
    TypeName("linear");


    // Constructors

        //- Construct from mesh
        linear(const fvMesh& mesh)
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}

        //- Construct from Istream
        linear
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream&
        )
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}

        //- Construct from faceFlux and Istream
        linear
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
        tmp<surfaceScalarField> weights(const VolField<Type>&) const
        {
            return this->mesh().surfaceInterpolation::weights();
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const linear&) = delete;
};


template<class Type>
tmp<SurfaceField<Type>>
linearInterpolate(const VolField<Type>& vf)
{
    return surfaceInterpolationScheme<Type>::interpolate
    (
        vf,
        vf.mesh().surfaceInterpolation::weights()
    );
}


template<class Type>
tmp<SurfaceField<Type>>
linearInterpolate(const tmp<VolField<Type>>& tvf)
{
    tmp<SurfaceField<Type>> tinterp = linearInterpolate(tvf());
    tvf.clear();
    return tinterp;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
