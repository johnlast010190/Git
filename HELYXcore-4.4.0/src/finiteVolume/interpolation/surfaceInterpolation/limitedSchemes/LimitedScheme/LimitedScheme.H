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
    (c) 2024 Engys Ltd.

Class
    Foam::LimitedScheme

Group
    grpFvLimitedSurfaceInterpolationSchemes

Description
    Class to create NVD/TVD limited weighting-factors.

    The particular differencing scheme class is supplied as a template
    argument, the weight function of which is called by the weight function
    of this class for the internal faces as well as faces of coupled
    patches (e.g. processor-processor patches). The weight function is
    supplied the central-differencing weighting factor, the face-flux, the
    cell and face gradients (from which the normalised variable
    distribution may be created) and the cell centre distance.

    This code organisation is both neat and efficient, allowing for
    convenient implementation of new schemes to run on parallelised cases.

SourceFiles
    LimitedScheme.C

\*---------------------------------------------------------------------------*/

#ifndef LimitedScheme_H
#define LimitedScheme_H

#include "interpolation/surfaceInterpolation/limitedSchemes/limitedSurfaceInterpolationScheme/limitedSurfaceInterpolationScheme.H"
#include "interpolation/surfaceInterpolation/limitedSchemes/LimitedScheme/LimitFuncs.H"
#include "interpolation/surfaceInterpolation/limitedSchemes/LimitedScheme/NVDTVD.H"
#include "interpolation/surfaceInterpolation/limitedSchemes/LimitedScheme/NVDVTVDV.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class LimitedScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class Limiter, template<class> class LimitFunc>
class LimitedScheme
:
    public limitedSurfaceInterpolationScheme<Type>,
    public Limiter
{
    // Private Member Functions

        //- Calculate the limiter
        void calcLimiter
        (
            const VolField<Type>& phi,
            surfaceScalarField& limiterField
        ) const;


public:

    //- Runtime type information
    TypeName("LimitedScheme");

    typedef Limiter LimiterType;

    // Constructors

        //- Construct from mesh and faceFlux and limiter scheme
        LimitedScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            const Limiter& weight
        )
        :
            limitedSurfaceInterpolationScheme<Type>(mesh, faceFlux),
            Limiter(weight)
        {}

        //- Construct from mesh and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        LimitedScheme
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& is
        )
        :
            limitedSurfaceInterpolationScheme<Type>(mesh, db, is),
            Limiter(is)
        {}

        //- Construct from mesh, faceFlux and Istream
        LimitedScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            limitedSurfaceInterpolationScheme<Type>(mesh, faceFlux),
            Limiter(is)
        {}

        //- Disallow default bitwise copy construction
        LimitedScheme(const LimitedScheme&) = delete;


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> limiter(const VolField<Type>&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LimitedScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeLimitedSurfaceInterpolationTypeScheme\
(                                                                              \
    SS,                                                                        \
    LIMITER,                                                                   \
    NVDTVD,                                                                    \
    LIMFUNC,                                                                   \
    TYPE                                                                       \
)                                                                              \
                                                                               \
typedef LimitedScheme<TYPE, LIMITER<NVDTVD>, limitFuncs::LIMFUNC>              \
    LimitedScheme##TYPE##LIMITER##NVDTVD##LIMFUNC##_;                          \
defineTemplateTypeNameAndDebugWithName                                         \
    (LimitedScheme##TYPE##LIMITER##NVDTVD##LIMFUNC##_, #SS, 0);                \
                                                                               \
surfaceInterpolationScheme<TYPE>::addMeshConstructorToTable                    \
<LimitedScheme<TYPE, LIMITER<NVDTVD>, limitFuncs::LIMFUNC>>                    \
    add##SS##LIMFUNC##TYPE##MeshConstructorToTable_;                           \
                                                                               \
surfaceInterpolationScheme<TYPE>::addMeshFluxConstructorToTable                \
<LimitedScheme<TYPE, LIMITER<NVDTVD>, limitFuncs::LIMFUNC>>                    \
    add##SS##LIMFUNC##TYPE##MeshFluxConstructorToTable_;                       \
                                                                               \
limitedSurfaceInterpolationScheme<TYPE>::addMeshConstructorToTable             \
<LimitedScheme<TYPE, LIMITER<NVDTVD>, limitFuncs::LIMFUNC>>                    \
    add##SS##LIMFUNC##TYPE##MeshConstructorToLimitedTable_;                    \
                                                                               \
limitedSurfaceInterpolationScheme<TYPE>::addMeshFluxConstructorToTable         \
<LimitedScheme<TYPE, LIMITER<NVDTVD>, limitFuncs::LIMFUNC>>                    \
    add##SS##LIMFUNC##TYPE##MeshFluxConstructorToLimitedTable_;


#define makeLimitedSurfaceInterpolationScheme(SS, LIMITER)                     \
                                                                               \
makeLimitedSurfaceInterpolationTypeScheme(SS,LIMITER,NVDTVD,magSqr,scalar)     \
makeLimitedSurfaceInterpolationTypeScheme(SS,LIMITER,NVDTVD,magSqr,vector)     \
makeLimitedSurfaceInterpolationTypeScheme                                      \
(                                                                              \
    SS,                                                                        \
    LIMITER,                                                                   \
    NVDTVD,                                                                    \
    magSqr,                                                                    \
    sphericalTensor                                                            \
)                                                                              \
makeLimitedSurfaceInterpolationTypeScheme(SS,LIMITER,NVDTVD,magSqr,symmTensor)\
makeLimitedSurfaceInterpolationTypeScheme(SS,LIMITER,NVDTVD,magSqr,tensor)


#define makeLimitedVSurfaceInterpolationScheme(SS, LIMITER)                    \
makeLimitedSurfaceInterpolationTypeScheme(SS,LIMITER,NVDVTVDV,null,vector)


#define makeLLimitedSurfaceInterpolationTypeScheme\
(                                                                              \
    SS,                                                                        \
    LLIMITER,                                                                  \
    LIMITER,                                                                   \
    NVDTVD,                                                                    \
    LIMFUNC,                                                                   \
    TYPE                                                                       \
)                                                                              \
                                                                               \
typedef LimitedScheme<TYPE, LLIMITER<LIMITER<NVDTVD>>, limitFuncs::LIMFUNC>    \
    LimitedScheme##TYPE##LLIMITER##LIMITER##NVDTVD##LIMFUNC##_;                \
defineTemplateTypeNameAndDebugWithName                                         \
    (LimitedScheme##TYPE##LLIMITER##LIMITER##NVDTVD##LIMFUNC##_, #SS, 0);      \
                                                                               \
surfaceInterpolationScheme<TYPE>::addMeshConstructorToTable                    \
<LimitedScheme<TYPE, LLIMITER<LIMITER<NVDTVD>>, limitFuncs::LIMFUNC>>          \
    add##SS##LIMFUNC##TYPE##MeshConstructorToTable_;                           \
                                                                               \
surfaceInterpolationScheme<TYPE>::addMeshFluxConstructorToTable                \
<LimitedScheme<TYPE, LLIMITER<LIMITER<NVDTVD>>, limitFuncs::LIMFUNC>>          \
    add##SS##LIMFUNC##TYPE##MeshFluxConstructorToTable_;                       \
                                                                               \
limitedSurfaceInterpolationScheme<TYPE>::addMeshConstructorToTable             \
<LimitedScheme<TYPE, LLIMITER<LIMITER<NVDTVD>>, limitFuncs::LIMFUNC>>          \
    add##SS##LIMFUNC##TYPE##MeshConstructorToLimitedTable_;                    \
                                                                               \
limitedSurfaceInterpolationScheme<TYPE>::addMeshFluxConstructorToTable         \
<LimitedScheme<TYPE, LLIMITER<LIMITER<NVDTVD>>, limitFuncs::LIMFUNC>>          \
    add##SS##LIMFUNC##TYPE##MeshFluxConstructorToLimitedTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interpolation/surfaceInterpolation/limitedSchemes/LimitedScheme/LimitedScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
