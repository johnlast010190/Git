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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2019-2024 Engys Ltd.

Class
    Foam::fv::gaussLaplacianScheme

Group
    grpFvLaplacianSchemes

Description
    Basic second-order laplacian using face-gradients and Gauss' theorem.

SourceFiles
    gaussLaplacianScheme.C

\*---------------------------------------------------------------------------*/

#ifndef gaussLaplacianScheme_H
#define gaussLaplacianScheme_H

#include "finiteVolume/laplacianSchemes/laplacianScheme/laplacianScheme.H"
#include "fvMatrices/fvBlockMatrix/fvBlockMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class gaussLaplacianScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GType>
class gaussLaplacianScheme
:
    public fv::laplacianScheme<Type, GType>
{
    // Private Member Functions

        tmp<SurfaceField<Type>> gammaSnGradCorr
        (
            const surfaceVectorField& SfGammaCorr,
            const VolField<Type>&
        );


public:

    //- Runtime type information
    TypeName("Gauss");


    // Constructors

        //- Construct null
        gaussLaplacianScheme(const fvMesh& mesh)
        :
            laplacianScheme<Type, GType>(mesh)
        {}

        //- Construct from Istream
        gaussLaplacianScheme
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& is,
            const word& forceGradSchemeName=word::null
        )
        :
            laplacianScheme<Type, GType>(mesh, db, is, forceGradSchemeName)
        {}

        //- Construct from mesh, interpolation and snGradScheme schemes
        gaussLaplacianScheme
        (
            const fvMesh& mesh,
            const tmp<surfaceInterpolationScheme<GType>>& igs,
            const tmp<snGradScheme<Type>>& sngs
        )
        :
            laplacianScheme<Type, GType>(mesh, igs, sngs)
        {}

        //- Disallow default bitwise copy construction
        gaussLaplacianScheme(const gaussLaplacianScheme&) = delete;


    //- Destructor
    virtual ~gaussLaplacianScheme()
    {}


    // Member Functions

        static tmp<fvMatrix<Type>> fvmLaplacianUncorrected
        (
            const surfaceScalarField& gammaMagSf,
            const surfaceScalarField& deltaCoeffs,
            const VolField<Type>&
        );

        tmp<VolField<Type>> fvcLaplacian
        (
            const VolField<Type>&
        );

        tmp<fvMatrix<Type>> fvmLaplacian
        (
            const SurfaceField<GType>&,
            const VolField<Type>&
        );

        tmp
        <
            fvBlockMatrix<Type>
        > fvmBLaplacian
        (
            const SurfaceField<GType>&,
            const VolField<Type>&
        );

        tmp<VolField<Type>> fvcLaplacian
        (
            const SurfaceField<GType>&,
            const VolField<Type>&
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const gaussLaplacianScheme&) = delete;
};




// Use macros to emulate partial-specialisation of the the Laplacian functions
// for scalar diffusivity gamma

#define defineFvmLaplacianScalarGamma(Type)                                    \
                                                                               \
template<>                                                                     \
inline tmp<fvMatrix<Type>> gaussLaplacianScheme<Type, scalar>::fvmLaplacian    \
(                                                                              \
    const SurfaceField<scalar>&,                 \
    const VolField<Type>&                         \
);                                                                             \
                                                                               \
template<>                                                                     \
tmp<VolField<Type>>                               \
inline gaussLaplacianScheme<Type, scalar>::fvcLaplacian                        \
(                                                                              \
    const SurfaceField<scalar>&,                 \
    const VolField<Type>&                         \
);


defineFvmLaplacianScalarGamma(scalar);
defineFvmLaplacianScalarGamma(vector);
defineFvmLaplacianScalarGamma(sphericalTensor);
defineFvmLaplacianScalarGamma(symmTensor);
defineFvmLaplacianScalarGamma(tensor);



#define defineFvmScalarBLaplacianGamma(Type)                                   \
                                                                               \
template<>                                                                     \
inline tmp<fvBlockMatrix<scalar>> gaussLaplacianScheme<scalar, Type>::fvmBLaplacian   \
(                                                                              \
    const SurfaceField<Type>&,                   \
    const VolField<scalar>&                       \
);                                                                             \


defineFvmScalarBLaplacianGamma(sphericalTensor);
defineFvmScalarBLaplacianGamma(symmTensor);
defineFvmScalarBLaplacianGamma(tensor);


#define defineFvmVectorBLaplacianGamma(Type)                                   \
                                                                               \
template<>                                                                     \
inline tmp<fvBlockMatrix<vector>> gaussLaplacianScheme<vector, Type>::fvmBLaplacian   \
(                                                                              \
    const SurfaceField<Type>&,                   \
    const VolField<vector>& vf                    \
)                                                                              \
{                                                                              \
    tmp<fvBlockMatrix<vector>> tfvbm                                           \
    (                                                                          \
        new fvBlockMatrix<vector>                                              \
        (                                                                      \
            const_cast<VolField<vector>&>(vf)     \
        )                                                                      \
    );                                                                         \
    FatalErrorInFunction                                                       \
        << " BLaplacianscalar(scalar, scalar-vector) suppoted"                 \
        << abort(FatalError);                                                  \
    return tfvbm;                                                              \
}                                                                              \


defineFvmVectorBLaplacianGamma(sphericalTensor);
defineFvmVectorBLaplacianGamma(symmTensor);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/laplacianSchemes/gaussLaplacianScheme/gaussLaplacianScheme.C"
    #include "finiteVolume/laplacianSchemes/gaussLaplacianScheme/gaussBlockLaplacianScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
