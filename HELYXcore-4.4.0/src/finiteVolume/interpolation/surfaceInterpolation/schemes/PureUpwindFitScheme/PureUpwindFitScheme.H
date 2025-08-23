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
    Foam::PureUpwindFitScheme

Group
    grpFvSurfaceInterpolationSchemes

Description
    Upwind biased fit surface interpolation scheme that applies an explicit
    correction to upwind.

\*---------------------------------------------------------------------------*/

#ifndef PureUpwindFitScheme_H
#define PureUpwindFitScheme_H

#include "interpolation/surfaceInterpolation/schemes/UpwindFitScheme/UpwindFitData.H"
#include "interpolation/surfaceInterpolation/limitedSchemes/upwind/upwind.H"
#include "primitives/bools/Switch/Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class PureUpwindFitScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class Polynomial, class Stencil>
class PureUpwindFitScheme
:
    public upwind<Type>
{
    // Private Data

        //- Factor the fit is allowed to deviate from linear.
        //  This limits the amount of high-order correction and increases
        //  stability on bad meshes
        const scalar linearLimitFactor_;

        //- Weights for central stencil
        const scalar centralWeight_;


public:

    //- Runtime type information
    TypeName("PureUpwindFitScheme");


    // Constructors

        //- Construct from mesh and Istream
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        PureUpwindFitScheme
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& is
        )
        :
            upwind<Type>
            (
                mesh,
                mesh.lookupObject<surfaceScalarField>(word(is))
            ),
            linearLimitFactor_(readScalar(is)),
            centralWeight_(1000)
        {}


        //- Construct from mesh, faceFlux and Istream
        PureUpwindFitScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            upwind<Type>(mesh, faceFlux),
            linearLimitFactor_(readScalar(is)),
            centralWeight_(1000)
        {}

        //- Disallow default bitwise copy construction
        PureUpwindFitScheme(const PureUpwindFitScheme&) = delete;


    // Member Functions

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<SurfaceField<Type>>
        correction(const VolField<Type>& vf) const
        {
            const fvMesh& mesh = this->mesh();

            // Use the owner/neighbour splitting constructor
            const extendedUpwindCellToFaceStencil& stencil = Stencil::New(mesh);

            // TODO: Make the `New` functions not be insane. Allocation routines should not
            //       return *references* omg.
#if defined(__clang__)
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdangling-reference"
#endif
            const UpwindFitData<Polynomial>& ufd =
            UpwindFitData<Polynomial>::New
            (
                mesh,
                stencil,
                false,              //offset to upwind
                linearLimitFactor_,
                centralWeight_
            );
#if defined(__clang__)
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

            const List<scalarList>& fo = ufd.owncoeffs();
            const List<scalarList>& fn = ufd.neicoeffs();

            return stencil.weightedSum(this->faceFlux_, vf, fo, fn);
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const PureUpwindFitScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makePureUpwindFitSurfaceInterpolationTypeScheme\
(                                                                              \
    SS,                                                                        \
    POLYNOMIAL,                                                                \
    STENCIL,                                                                   \
    TYPE                                                                       \
)                                                                              \
                                                                               \
typedef PureUpwindFitScheme<TYPE, POLYNOMIAL, STENCIL>                         \
    PureUpwindFitScheme##TYPE##POLYNOMIAL##STENCIL##_;                         \
defineTemplateTypeNameAndDebugWithName                                         \
    (PureUpwindFitScheme##TYPE##POLYNOMIAL##STENCIL##_, #SS, 0);               \
                                                                               \
surfaceInterpolationScheme<TYPE>::addMeshConstructorToTable                    \
<PureUpwindFitScheme<TYPE, POLYNOMIAL, STENCIL>>                               \
    add##SS##STENCIL##TYPE##MeshConstructorToTable_;                           \
                                                                               \
surfaceInterpolationScheme<TYPE>::addMeshFluxConstructorToTable                \
<PureUpwindFitScheme<TYPE, POLYNOMIAL, STENCIL>>                               \
    add##SS##STENCIL##TYPE##MeshFluxConstructorToTable_;

#define makePureUpwindFitSurfaceInterpolationScheme(SS, POLYNOMIAL, STENCIL)   \
                                                                               \
makePureUpwindFitSurfaceInterpolationTypeScheme(SS,POLYNOMIAL,STENCIL,scalar) \
makePureUpwindFitSurfaceInterpolationTypeScheme(SS,POLYNOMIAL,STENCIL,vector) \
makePureUpwindFitSurfaceInterpolationTypeScheme                                \
(                                                                              \
    SS,                                                                        \
    POLYNOMIAL,                                                                \
    STENCIL,                                                                   \
    sphericalTensor                                                            \
)                                                                              \
makePureUpwindFitSurfaceInterpolationTypeScheme                                \
(                                                                              \
    SS,                                                                        \
    POLYNOMIAL,                                                                \
    STENCIL,                                                                   \
    symmTensor                                                                 \
)                                                                              \
makePureUpwindFitSurfaceInterpolationTypeScheme(SS,POLYNOMIAL,STENCIL,tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
