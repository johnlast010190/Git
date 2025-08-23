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
    Foam::surfaceInterpolationScheme

Description
    Abstract base class for surface interpolation schemes.

SourceFiles
    surfaceInterpolationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceInterpolationScheme_H
#define surfaceInterpolationScheme_H

#include "memory/tmp/tmp.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "db/typeInfo/typeInfo.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "db/objectRegistry/objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                 Class surfaceInterpolationScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class surfaceInterpolationScheme
:
    public tmp<surfaceInterpolationScheme<Type>>::refCount
{
    // Private data

        //- Hold reference to mesh
        const fvMesh& mesh_;


public:

    //- Runtime type information
    TypeName("surfaceInterpolationScheme");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            surfaceInterpolationScheme,
            Mesh,
            (
                const fvMesh& mesh,
                const objectRegistry& db,
                Istream& schemeData
            ),
            (mesh, db, schemeData)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            surfaceInterpolationScheme,
            MeshFlux,
            (
                const fvMesh& mesh,
                const surfaceScalarField& faceFlux,
                Istream& schemeData
            ),
            (mesh, faceFlux, schemeData)
        );


    // Constructors

        //- Construct from mesh
        surfaceInterpolationScheme(const fvMesh& mesh)
        :
            mesh_(mesh)
        {}

        //- Disallow default bitwise copy construction
        surfaceInterpolationScheme(const surfaceInterpolationScheme&) = delete;


    // Selectors

        //- Return new tmp interpolation scheme
        static tmp<surfaceInterpolationScheme<Type>> New
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& schemeData
        );

        //- Return new tmp interpolation scheme
        static tmp<surfaceInterpolationScheme<Type>> New
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        );


    //- Destructor
    virtual ~surfaceInterpolationScheme();


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }


        //- Return the face-interpolate of the given cell field
        //  with the given owner and neighbour weigting factors
        static tmp<SurfaceField<Type>>
        interpolate
        (
            const VolField<Type>&,
            const tmp<surfaceScalarField>&,
            const tmp<surfaceScalarField>&
        );

        //- Return the face-interpolate of the given cell field
        //  with the given weighting factors dotted with given field Sf
        template<class SFType>
        static tmp
        <
            SurfaceField
            <
                typename innerProduct<typename SFType::value_type, Type>::type
            >
        >
        dotInterpolate
        (
            const SFType& Sf,
            const VolField<Type>& vf,
            const tmp<surfaceScalarField>& tlambdas
        );

        //- Return the face-interpolate of the given cell field
        //  with the given weighting factors
        static tmp<SurfaceField<Type>>
        interpolate
        (
            const VolField<Type>&,
            const tmp<surfaceScalarField>&
        );

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> weights
        (
            const VolField<Type>&
        ) const = 0;

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return false;
        }

        //- Return the explicit correction to the face-interpolate
        //  for the given field
        virtual tmp<SurfaceField<Type>> correction(const VolField<Type>&) const
        {
            return tmp<SurfaceField<Type>>(nullptr);
        }

        //- Return the face-interpolate of the given cell field
        //  with explicit correction dotted with given field Sf
        virtual tmp<SurfaceField<typename innerProduct<vector, Type>::type>>
        dotInterpolate
        (
            const surfaceVectorField& Sf,
            const VolField<Type>& vf
        ) const;

        //- Return the face-interpolate of the given tmp cell field
        //  with explicit correction dotted with given field Sf
        tmp<SurfaceField<typename innerProduct<vector, Type>::type>>
        dotInterpolate
        (
            const surfaceVectorField& Sf,
            const tmp<VolField<Type>>&
        ) const;

        //- Return the face-interpolate of the given cell field
        //  with explicit correction
        virtual tmp<SurfaceField<Type>>
        interpolate(const VolField<Type>&) const;

        //- Return the face-interpolate of the given tmp cell field
        //  with explicit correction
        tmp<SurfaceField<Type>> interpolate(const tmp<VolField<Type>>&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const surfaceInterpolationScheme&) = delete;
};


template<>
tmp
<
    GeometricField
    <
        typename innerProduct<vector, scalar>::type,
        fvsPatchField,
        surfaceMesh
    >
>
surfaceInterpolationScheme<scalar>::dotInterpolate
(
    const surfaceVectorField& Sf,
    const VolField<scalar>&
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeSurfaceInterpolationTypeScheme(SS, Type)                           \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
surfaceInterpolationScheme<Type>::addMeshConstructorToTable<SS<Type>>          \
    add##SS##Type##MeshConstructorToTable_;                                    \
                                                                               \
surfaceInterpolationScheme<Type>::addMeshFluxConstructorToTable<SS<Type>>      \
    add##SS##Type##MeshFluxConstructorToTable_;

#define makeSurfaceInterpolationScheme(SS)                                     \
                                                                               \
makeSurfaceInterpolationTypeScheme(SS, scalar)                                 \
makeSurfaceInterpolationTypeScheme(SS, vector)                                 \
makeSurfaceInterpolationTypeScheme(SS, sphericalTensor)                        \
makeSurfaceInterpolationTypeScheme(SS, symmTensor)                             \
makeSurfaceInterpolationTypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
