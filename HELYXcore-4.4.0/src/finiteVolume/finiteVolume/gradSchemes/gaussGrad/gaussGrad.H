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
    (c) 2017-2019 Engys Ltd.

Class
    Foam::fv::gaussGrad

Group
    grpFvGradSchemes

Description
    Basic second-order gradient scheme using face-interpolation
    and Gauss' theorem.

SourceFiles
    gaussGrad.C

\*---------------------------------------------------------------------------*/

#ifndef gaussGrad_H
#define gaussGrad_H

#include "finiteVolume/gradSchemes/gradScheme/gradScheme.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.H"
#include "interpolation/surfaceInterpolation/schemes/linear/linear.H"
#include "blockLduSystem/BlockLduSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class gaussGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class gaussGrad
:
    public fv::gradScheme<Type>
{

protected:

    // Protected Data

        tmp<surfaceInterpolationScheme<Type>> tinterpScheme_;


public:

    //- Runtime type information
    TypeName("Gauss");


    // Constructors

        //- Construct from mesh
        gaussGrad(const fvMesh& mesh)
        :
            gradScheme<Type>(mesh),
            tinterpScheme_(new linear<Type>(mesh))
        {}

        //- Construct from mesh and Istream
        gaussGrad
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& is
        )
        :
            gradScheme<Type>(mesh),
            tinterpScheme_(nullptr)
        {
            if (is.eof())
            {
                tinterpScheme_ =
                    tmp<surfaceInterpolationScheme<Type>>
                    (
                        new linear<Type>(mesh)
                    );
            }
            else
            {
                tinterpScheme_ =
                    tmp<surfaceInterpolationScheme<Type>>
                    (
                        surfaceInterpolationScheme<Type>::New(mesh, db, is)
                    );
            }
        }

        //- Disallow default bitwise copy construction
        gaussGrad(const gaussGrad&) = delete;


    // Member Functions

        //- Surface field internal face gradient contribution
        static tmp<VolField<typename outerProduct<vector, Type>::type>>
        gradFinternal
        (
            const SurfaceField<Type>&,
            const word& name
        );

        //- Return the gradient of the given field
        //  calculated using Gauss' theorem on the given surface field
        static tmp<VolField<typename outerProduct<vector, Type>::type>> gradf
        (
            const SurfaceField<Type>&,
            const word& name
        );

        //- Return the gradient of the given field to the gradScheme::grad
        //  for optional caching
        virtual tmp<VolField<typename outerProduct<vector, Type>::type>>
        calcGrad
        (
            const VolField<Type>& vsf,
            const word& name
        ) const;

        tmp
        <
            BlockLduSystem<vector, typename outerProduct<vector, Type>::type>
        > fvmGrad
        (
            const VolField<Type>&
        ) const;

        //- Correct the boundary values of the gradient using the patchField
        // snGrad functions
        static void correctBoundaryConditions
        (
            const VolField<Type>&,
            VolField<typename outerProduct<vector, Type>::type>&
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const gaussGrad&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/gradSchemes/gaussGrad/gaussGrad.C"
#endif

// Specialisation for scalars
#include "finiteVolume/gradSchemes/gaussGrad/scalarGaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
