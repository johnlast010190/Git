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
    Foam::clippedReverseLinear

Group
    grpFvSurfaceInterpolationSchemes

Description
    Inversed weight central-differencing interpolation scheme class
    with clipping/limiter.

SourceFiles
    clippedReverseLinear.C

\*---------------------------------------------------------------------------*/

#ifndef clippedReverseLinear_H
#define clippedReverseLinear_H

#include "interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class clippedReverseLinear Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class clippedReverseLinear
:
    public surfaceInterpolationScheme<Type>
{
    // Private data

        const scalar cellSizeRatio_;
        scalar wfLimit_;


    // Private Member Functions

        void calcWfLimit()
        {
            if (cellSizeRatio_ <= 0 || cellSizeRatio_ > 1)
            {
                FatalErrorInFunction
                    << "Given cellSizeRatio of " << cellSizeRatio_
                    << " is not between 0 and 1"
                    << exit(FatalError);
            }

            wfLimit_ = cellSizeRatio_/(1.0 + cellSizeRatio_);
        }


public:

    //- Runtime type information
    TypeName("clippedReverseLinear");


    // Constructors

        //- Construct from mesh and cellSizeRatio
        clippedReverseLinear(const fvMesh& mesh, const scalar cellSizeRatio)
        :
            surfaceInterpolationScheme<Type>(mesh),
            cellSizeRatio_(cellSizeRatio)
        {
            calcWfLimit();
        }

        //- Construct from Istream
        clippedReverseLinear
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            cellSizeRatio_(readScalar(is))
        {
            calcWfLimit();
        }

        //- Construct from faceFlux and Istream
        clippedReverseLinear
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            cellSizeRatio_(readScalar(is))
        {
            calcWfLimit();
        }


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const VolField<Type>&
        ) const
        {
            const fvMesh& mesh = this->mesh();

            tmp<surfaceScalarField> tcdWeights
            (
                mesh.surfaceInterpolation::weights()
            );
            const surfaceScalarField& cdWeights = tcdWeights();

            tmp<surfaceScalarField> tclippedReverseLinearWeights
            (
                surfaceScalarField::New
                (
                    "clippedReverseLinearWeights",
                    1.0 - cdWeights
                )
            );
            surfaceScalarField& clippedReverseLinearWeights =
                tclippedReverseLinearWeights.ref();

            const scalar& minWei = wfLimit_;
            const scalar maxWei = 1.0 - minWei;

            clippedReverseLinearWeights.primitiveFieldRef() =
                max
                (
                    min
                    (
                        clippedReverseLinearWeights.primitiveField(),
                        maxWei
                    ),
                    minWei
                );

            surfaceScalarField::Boundary& clwbf =
                clippedReverseLinearWeights.boundaryFieldRef();

            forAll(mesh.boundary(), patchi)
            {
                if (clwbf[patchi].coupled())
                {
                    clwbf[patchi] =
                        max
                        (
                            min
                            (
                                clwbf[patchi],
                                maxWei
                            ),
                            minWei
                        );
                }
                else
                {
                    clwbf[patchi] = cdWeights.boundaryField()[patchi];
                }
            }

            return tclippedReverseLinearWeights;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const clippedReverseLinear&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
