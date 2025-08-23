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
    Foam::multivariateSurfaceInterpolationScheme

Description
    Abstract base class for multi-variate surface interpolation schemes.

SourceFiles
    multivariateSurfaceInterpolationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef multivariateSurfaceInterpolationScheme_H
#define multivariateSurfaceInterpolationScheme_H

#include "interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.H"
#include "containers/HashTables/HashTable/HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class multivariateSurfaceInterpolationScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class multivariateSurfaceInterpolationScheme
:
    public tmp<multivariateSurfaceInterpolationScheme<Type>>::refCount
{

public:

    //- fieldTable
    class fieldTable
    :
        public HashTable<const VolField<Type>*>
    {
    public:

        fieldTable()
        {}

        void add(const VolField<Type>& f)
        {
            this->insert(f.name(), &f);
        }
    };


private:

    // Private data

        //- Hold reference to mesh
        const fvMesh& mesh_;

        //- HashTable of pointers to the field set
        const fieldTable& fields_;


    // Private Member Functions


public:

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            multivariateSurfaceInterpolationScheme,
            Istream,
            (
                const fvMesh& mesh,
                const fieldTable& fields,
                const surfaceScalarField& faceFlux,
                Istream& is
            ),
            (mesh, fields, faceFlux, is)
        );


    // Constructors

        //- Construct for interpolating given field
        multivariateSurfaceInterpolationScheme
        (
            const fvMesh& mesh,
            const fieldTable& fields,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        );

        //- Disallow default bitwise copy construction
        multivariateSurfaceInterpolationScheme
        (
            const multivariateSurfaceInterpolationScheme&
        ) = delete;


    // Selectors

        //- Return a pointer to a new gradScheme created on freestore
        static tmp<multivariateSurfaceInterpolationScheme<Type>> New
        (
            const fvMesh& mesh,
            const fieldTable& fields,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        );


    //- Destructor
    virtual ~multivariateSurfaceInterpolationScheme();


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return fields to be interpolated
        const fieldTable& fields() const
        {
            return fields_;
        }


    // Member Operators

        //- surfaceInterpolationScheme sub-class returned by operator(field)
        //  which is used as the interpolation scheme for the field
        class fieldScheme
        :
            public surfaceInterpolationScheme<Type>
        {

        public:

            // Constructors

                fieldScheme
                (
                    const VolField<Type>& field
                )
                :
                    surfaceInterpolationScheme<Type>(field.mesh())
                {}


            // Member Functions

                //- Return the interpolation weighting factors
                virtual tmp<surfaceScalarField> weights
                (
                    const VolField<Type>& field
                ) const = 0;
        };

        virtual tmp<surfaceInterpolationScheme<Type>> operator()
        (
            const VolField<Type>& field
        ) const = 0;

        //- Disallow default bitwise assignment
        void operator=(const multivariateSurfaceInterpolationScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeMultivariateSurfaceInterpolationTypeScheme(SS, Type)               \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
multivariateSurfaceInterpolationScheme<Type>::                                 \
addIstreamConstructorToTable<SS<Type>>                                         \
    add##SS##Type##ConstructorToTable_;


#define makeMultivariateSurfaceInterpolationScheme(SS)                         \
                                                                               \
makeMultivariateSurfaceInterpolationTypeScheme(SS, scalar)                     \
makeMultivariateSurfaceInterpolationTypeScheme(SS, vector)                     \
makeMultivariateSurfaceInterpolationTypeScheme(SS, sphericalTensor)            \
makeMultivariateSurfaceInterpolationTypeScheme(SS, symmTensor)                 \
makeMultivariateSurfaceInterpolationTypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interpolation/surfaceInterpolation/multivariateSchemes/multivariateSurfaceInterpolationScheme/multivariateSurfaceInterpolationScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
