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
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2024-2025 Engys Ltd.

Class
    Foam::fv::localEulerDdtScheme

Group
    grpFvDdtSchemes

Description
    Local time-step first-order Euler implicit/explicit ddt.

    The reciprocal of the local time-step field is looked-up from the database.

    This scheme should only be used for steady-state computations using
    transient codes where local time-stepping is preferably to under-relaxation
    for transport consistency reasons.

See also
    Foam::fv::CoEulerDdtScheme

SourceFiles
    localEulerDdt.C
    localEulerDdtScheme.C
    localEulerDdtSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef localEulerDdtScheme_H
#define localEulerDdtScheme_H

#include "finiteVolume/ddtSchemes/ddtScheme/ddtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class localEulerDdt Declaration
\*---------------------------------------------------------------------------*/

class localEulerDdt
{
public:

    //- Name of the reciprocal local time-step field
    static word rDeltaTName;

    //- Name of the reciprocal local face time-step field
    static word rDeltaTfName;

    //- Name of the reciprocal local sub-cycling time-step field
    static word rSubDeltaTName;


    // Constructors

        localEulerDdt()
        {}


    // Member Functions

        //- Return true if LTS is enabled
        static bool enabled(const fvMesh& mesh);

        //- Return the reciprocal of the local time-step
        //  looked-up from the objectRegistry
        static const volScalarField& localRDeltaT(const fvMesh& mesh);

        //- Return the reciprocal of the local face time-step
        //  looked-up from the objectRegistry
        static const surfaceScalarField& localRDeltaTf(const fvMesh& mesh);

        //- Calculate and return the reciprocal of the local sub-cycling
        //  time-step
        static tmp<volScalarField> localRSubDeltaT
        (
            const fvMesh& mesh,
            const label nAlphaSubCycles
        );
};


/*---------------------------------------------------------------------------*\
                       Class localEulerDdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class localEulerDdtScheme
:
    public localEulerDdt,
    public fv::ddtScheme<Type>
{
    // Private Member Functions

        //- Return the reciprocal of the local time-step
        const volScalarField& localRDeltaT() const;

        //- Return the reciprocal of the local face time-step
        const surfaceScalarField& localRDeltaTf() const;


public:

    //- Runtime type information
    TypeName("localEuler");


    // Constructors

        //- Construct from mesh
        localEulerDdtScheme(const fvMesh& mesh)
        :
            ddtScheme<Type>(mesh)
        {}

        //- Construct from mesh and Istream
        localEulerDdtScheme(const fvMesh& mesh, Istream& is)
        :
            ddtScheme<Type>(mesh, is)
        {}

        //- Disallow default bitwise copy construction
        localEulerDdtScheme(const localEulerDdtScheme&) = delete;


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return fv::ddtScheme<Type>::mesh();
        }

        tmp<VolField<Type>> fvcDdt
        (
            const dimensioned<Type>&
        );

        tmp<VolField<Type>> fvcDdt
        (
            const VolField<Type>&
        );

        tmp<VolField<Type>> fvcDdt
        (
            const dimensionedScalar&,
            const VolField<Type>&
        );

        tmp<VolField<Type>> fvcDdt
        (
            const volScalarField&,
            const VolField<Type>&
        );

        tmp<VolField<Type>> fvcDdt
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const VolField<Type>& psi
        );

        tmp<SurfaceField<Type>> fvcDdt
        (
            const SurfaceField<Type>&
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const VolField<Type>&
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const dimensionedScalar&,
            const VolField<Type>&
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const volScalarField&,
            const VolField<Type>&
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const VolField<Type>& psi
        );

        typedef typename ddtScheme<Type>::fluxFieldType fluxFieldType;

        tmp<fluxFieldType> fvcDdtUfCorr
        (
            const VolField<Type>& U,
            const SurfaceField<Type>& Uf,
            const word
        );

        tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const VolField<Type>& U,
            const fluxFieldType& phi,
            const word
        );

        tmp<fluxFieldType> fvcDdtUfCorr
        (
            const volScalarField& rho,
            const VolField<Type>& U,
            const SurfaceField<Type>& Uf,
            const word
        );

        tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const volScalarField& rho,
            const VolField<Type>& U,
            const fluxFieldType& phi,
            const word
        );

        tmp<surfaceScalarField> meshPhi
        (
            const VolField<Type>&
        );

        tmp<scalarField> meshPhi(const VolField<Type>&, const label patchi);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const localEulerDdtScheme&) = delete;
};


template<>
tmp<surfaceScalarField> localEulerDdtScheme<scalar>::fvcDdtUfCorr
(
    const VolField<scalar>& U,
    const SurfaceField<scalar>& Uf,
    const word
);

template<>
tmp<surfaceScalarField> localEulerDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& U,
    const surfaceScalarField& phi,
    const word
);

template<>
tmp<surfaceScalarField> localEulerDdtScheme<scalar>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& Uf,
    const word
);

template<>
tmp<surfaceScalarField> localEulerDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& phi,
    const word
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/ddtSchemes/localEulerDdtScheme/localEulerDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
