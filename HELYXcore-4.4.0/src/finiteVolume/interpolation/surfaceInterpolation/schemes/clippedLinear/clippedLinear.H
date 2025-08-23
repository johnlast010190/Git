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
    Foam::clippedLinear

Group
    grpFvSurfaceInterpolationSchemes

Description
    Central-differencing interpolation scheme using clipped-weights to
    improve stability on meshes with very rapid variations in cell size.

SourceFiles
    clippedLinear.C

\*---------------------------------------------------------------------------*/

#ifndef clippedLinear_H
#define clippedLinear_H

#include "interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class clippedLinear Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class clippedLinear
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
    TypeName("clippedLinear");


    // Constructors

        //- Construct from mesh and cellSizeRatio
        clippedLinear(const fvMesh& mesh, const scalar cellSizeRatio)
        :
            surfaceInterpolationScheme<Type>(mesh),
            cellSizeRatio_(cellSizeRatio)
        {
            calcWfLimit();
        }

        //- Construct from Istream
        clippedLinear
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
        clippedLinear
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

            tmp<surfaceScalarField> tclippedLinearWeights
            (
                surfaceScalarField::New
                (
                    "clippedLinearWeights",
                    mesh,
                    dimless
                )
            );
            surfaceScalarField& clippedLinearWeights =
                tclippedLinearWeights.ref();

            clippedLinearWeights.primitiveFieldRef() =
                max(min(cdWeights.primitiveField(), 1 - wfLimit_), wfLimit_);

            surfaceScalarField::Boundary& clwbf =
                clippedLinearWeights.boundaryFieldRef();

            forAll(mesh.boundary(), patchi)
            {
                if (clwbf[patchi].coupled())
                {
                    clwbf[patchi] =
                        max
                        (
                            min
                            (
                                cdWeights.boundaryField()[patchi],
                                1 - wfLimit_
                            ),
                            wfLimit_
                        );
                }
                else
                {
                    clwbf[patchi] = cdWeights.boundaryField()[patchi];
                }
            }

            return tclippedLinearWeights;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const clippedLinear&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
