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
    (c) 2017 OpenCFD Ltd.
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2024-2025 Engys Ltd.

Class
    Foam::fv::ddtScheme

Group
    grpFvDdtSchemes

Description
    Abstract base class for ddt schemes.

SourceFiles
    ddtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef ddtScheme_H
#define ddtScheme_H

#include "memory/tmp/tmp.H"
#include "dimensionedTypes/dimensionedType/dimensionedType.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "db/typeInfo/typeInfo.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class fvMatrix;

class fvMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class ddtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class ddtScheme
:
    public tmp<ddtScheme<Type>>::refCount
{
protected:

    // Protected data

        const fvMesh& mesh_;

        //- Input for fvcDdtPhiCoeff (-1 default)
        autoPtr<Function1<scalar>> ddtPhiCoeff_;


    // Protected member function
        scalar getDdtPhiCoeff(const fvMesh& mesh);


public:

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            ddtScheme,
            Istream,
            (const fvMesh& mesh, Istream& schemeData),
            (mesh, schemeData)
        );


    // Constructors

        //- Construct from mesh
        ddtScheme(const fvMesh& mesh)
        :
            mesh_(mesh),
            ddtPhiCoeff_()
        {}

        //- Construct from mesh and Istream
        ddtScheme(const fvMesh& mesh, Istream& is)
        :
            mesh_(mesh),
            ddtPhiCoeff_()
        {}

        //- Disallow default bitwise copy construction
        ddtScheme(const ddtScheme&) = delete;


    // Selectors

        //- Return a pointer to a new ddtScheme created on freestore
        static tmp<ddtScheme<Type>> New
        (
            const fvMesh& mesh,
            Istream& schemeData
        );


    //- Destructor
    virtual ~ddtScheme();


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        virtual scalar pseudoDeltaTValue() const
        {
            return mesh().time().deltaTValue();
        }

        virtual tmp<VolField<Type>> fvcDdt
        (
            const dimensioned<Type>&
        ) = 0;

        virtual tmp<VolField<Type>> fvcDdt
        (
            const VolField<Type>&
        ) = 0;

        virtual tmp<VolField<Type>> fvcDdt
        (
            const dimensionedScalar&,
            const VolField<Type>&
        ) = 0;

        virtual tmp<VolField<Type>> fvcDdt
        (
            const volScalarField&,
            const VolField<Type>&
        ) = 0;

        virtual tmp<VolField<Type>> fvcDdt
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const VolField<Type>&
        ) = 0;

        virtual tmp<SurfaceField<Type>> fvcDdt
        (
            const SurfaceField<Type>&
        );

        virtual tmp<fvMatrix<Type>> fvmDdt
        (
            const VolField<Type>&
        ) = 0;

        virtual tmp<fvMatrix<Type>> fvmDdt
        (
            const dimensionedScalar&,
            const VolField<Type>&
        ) = 0;

        virtual tmp<fvMatrix<Type>> fvmDdt
        (
            const volScalarField&,
            const VolField<Type>&
        ) = 0;

        virtual tmp<fvMatrix<Type>> fvmDdt
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const VolField<Type>& vf
        ) = 0;

        typedef GeometricField
        <
            typename flux<Type>::type,
            fvsPatchField,
            surfaceMesh
        > fluxFieldType;

        tmp<surfaceScalarField> fvcDdtPhiCoeff
        (
            const VolField<Type>& U,
            const fluxFieldType& phi,
            const fluxFieldType& phiCorr
        );

        tmp<surfaceScalarField> fvcDdtPhiCoeff
        (
            const VolField<Type>& U,
            const fluxFieldType& phi
        );

        virtual tmp<fluxFieldType> fvcDdtUfCorr
        (
            const VolField<Type>& U,
            const SurfaceField<Type>& Uf,
            const word
        ) = 0;

        virtual tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const VolField<Type>& U,
            const fluxFieldType& phi,
            const word
        ) = 0;

        virtual tmp<fluxFieldType> fvcDdtUfCorr
        (
            const volScalarField& rho,
            const VolField<Type>& U,
            const SurfaceField<Type>& Uf,
            const word
        ) = 0;

        virtual tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const volScalarField& rho,
            const VolField<Type>& U,
            const fluxFieldType& phi,
            const word
        ) = 0;

        virtual tmp<surfaceScalarField> meshPhi
        (
            const VolField<Type>&
        ) = 0;

        virtual tmp<scalarField> meshPhi
        (
            const VolField<Type>&,
            const label patchi
        ) = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ddtScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<ddtScheme<Type>> getDdtScheme
(
    const VolField<Type>& f
)
{
    return
        ddtScheme<Type>::New
        (
            f.mesh(), f.mesh().schemes().ddtScheme("ddt(" + f.name() + ')')
        );
}

template<class FieldType, class Type>
static tmp<ddtScheme<Type>> getDdtScheme
(
    const FieldType& f1,
    const VolField<Type>& f
)
{
    return
        ddtScheme<Type>::New
        (
            f.mesh(),
            f.mesh().schemes().ddtScheme
            (
                "ddt("+f1.name()+','+f.name()+')'
            )
        );
}

template<class FieldType1, class FieldType2, class Type>
tmp<ddtScheme<Type>> getDdtScheme
(
    const FieldType1& f1,
    const FieldType2& f2,
    const VolField<Type>& f
)
{
    return
        ddtScheme<Type>::New
        (
            f.mesh(),
            f.mesh().schemes().ddtScheme
            (
                "ddt("+f1.name()+','+f2.name()+','+f.name()+')'
            )
        );
}

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeFvDdtTypeScheme(SS, Type)                                          \
    defineNamedTemplateTypeNameAndDebug(Foam::fv::SS<Foam::Type>, 0);          \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace fv                                                           \
        {                                                                      \
            ddtScheme<Type>::addIstreamConstructorToTable<SS<Type>>            \
                add##SS##Type##IstreamConstructorToTable_;                     \
        }                                                                      \
    }

#define makeFvDdtScheme(SS)                                                    \
                                                                               \
makeFvDdtTypeScheme(SS, scalar)                                                \
makeFvDdtTypeScheme(SS, vector)                                                \
makeFvDdtTypeScheme(SS, sphericalTensor)                                       \
makeFvDdtTypeScheme(SS, symmTensor)                                            \
makeFvDdtTypeScheme(SS, tensor)                                                \
                                                                               \
namespace Foam                                                                 \
{                                                                              \
namespace fv                                                                   \
{                                                                              \
                                                                               \
template<>                                                                     \
tmp<surfaceScalarField> SS<scalar>::fvcDdtUfCorr                               \
(                                                                              \
    const volScalarField& U,                                                   \
    const surfaceScalarField& Uf,                                              \
    const word interName                                                       \
)                                                                              \
{                                                                              \
    NotImplemented;                                                            \
}                                                                              \
                                                                               \
template<>                                                                     \
tmp<surfaceScalarField> SS<scalar>::fvcDdtPhiCorr                              \
(                                                                              \
    const volScalarField& U,                                                   \
    const surfaceScalarField& phi,                                             \
    const word interName                                                       \
)                                                                              \
{                                                                              \
    NotImplemented;                                                            \
}                                                                              \
                                                                               \
template<>                                                                     \
tmp<surfaceScalarField> SS<scalar>::fvcDdtUfCorr                               \
(                                                                              \
    const volScalarField& rho,                                                 \
    const volScalarField& U,                                                   \
    const surfaceScalarField& Uf,                                              \
    const word interName                                                       \
)                                                                              \
{                                                                              \
    NotImplemented;                                                            \
}                                                                              \
                                                                               \
template<>                                                                     \
tmp<surfaceScalarField> SS<scalar>::fvcDdtPhiCorr                              \
(                                                                              \
    const volScalarField& rho,                                                 \
    const volScalarField& U,                                                   \
    const surfaceScalarField& phi,                                             \
    const word interName                                                       \
)                                                                              \
{                                                                              \
    NotImplemented;                                                            \
}                                                                              \
                                                                               \
}                                                                              \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/ddtSchemes/ddtScheme/ddtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
