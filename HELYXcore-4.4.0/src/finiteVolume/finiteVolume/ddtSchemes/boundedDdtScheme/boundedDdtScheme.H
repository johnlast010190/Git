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
    (c) 2012-2022 OpenFOAM Foundation
    (c) 2024-2025 Engys Ltd.

Class
    Foam::fv::boundedDdtScheme

Group
    grpFvDdtSchemes

Description
    Bounded form of the selected ddt scheme.

    Boundedness is achieved by subtracting ddt(phi)*vf or Sp(ddt(rho), vf)
    which is non-conservative if ddt(rho) != 0 but conservative otherwise.

    Can be used for the ddt of bounded scalar properties to improve stability
    if insufficient convergence of the pressure equation causes temporary
    divergence of the flux field.

SourceFiles
    boundedDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef boundedDdtScheme_H
#define boundedDdtScheme_H

#include "finiteVolume/ddtSchemes/ddtScheme/ddtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class boundedDdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class boundedDdtScheme
:
    public fv::ddtScheme<Type>
{
    // Private data

        tmp<fv::ddtScheme<Type>> scheme_;


public:

    //- Runtime type information
    TypeName("bounded");


    // Constructors

        //- Construct from mesh and Istream
        boundedDdtScheme(const fvMesh& mesh, Istream& is)
        :
            ddtScheme<Type>(mesh, is),
            scheme_
            (
                fv::ddtScheme<Type>::New(mesh, is)
            )
        {}

        //- Disallow default bitwise copy construction
        boundedDdtScheme(const boundedDdtScheme&) = delete;


    // Member Functions

        //- Return mesh reference
        virtual const fvMesh& mesh() const
        {
            return fv::ddtScheme<Type>::mesh();
        }

        virtual scalar pseudoDeltaTValue() const
        {
            return scheme_().pseudoDeltaTValue();
        }

        virtual tmp<fv::ddtScheme<Type>> scheme() const
        {
            return tmp<fv::ddtScheme<Type>>(scheme_);
        }

        virtual tmp<VolField<Type>> fvcDdt
        (
            const dimensioned<Type>&
        );

        virtual tmp<VolField<Type>> fvcDdt
        (
            const VolField<Type>&
        );

        virtual tmp<VolField<Type>> fvcDdt
        (
            const dimensionedScalar&,
            const VolField<Type>&
        );

        virtual tmp<VolField<Type>> fvcDdt
        (
            const volScalarField&,
            const VolField<Type>&
        );

        virtual tmp<VolField<Type>> fvcDdt
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const VolField<Type>& psi
        );

        virtual tmp<fvMatrix<Type>> fvmDdt
        (
            const VolField<Type>&
        );

        virtual tmp<fvMatrix<Type>> fvmDdt
        (
            const dimensionedScalar&,
            const VolField<Type>&
        );

        virtual tmp<fvMatrix<Type>> fvmDdt
        (
            const volScalarField&,
            const VolField<Type>&
        );

        virtual tmp<fvMatrix<Type>> fvmDdt
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const VolField<Type>& psi
        );

        typedef typename ddtScheme<Type>::fluxFieldType fluxFieldType;

        virtual tmp<fluxFieldType> fvcDdtUfCorr
        (
            const VolField<Type>& U,
            const SurfaceField<Type>& Uf,
            const word
        );

        virtual tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const VolField<Type>& U,
            const fluxFieldType& phi,
            const word
        );

        virtual tmp<fluxFieldType> fvcDdtUfCorr
        (
            const volScalarField& rho,
            const VolField<Type>& U,
            const SurfaceField<Type>& Uf,
            const word
        );

        virtual tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const volScalarField& rho,
            const VolField<Type>& U,
            const fluxFieldType& phi,
            const word
        );

        virtual tmp<surfaceScalarField> meshPhi
        (
            const VolField<Type>&
        );

        virtual tmp<scalarField> meshPhi
        (
            const VolField<Type>&,
            const label patchi
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const boundedDdtScheme&) = delete;
};


template<>
tmp<surfaceScalarField> boundedDdtScheme<scalar>::fvcDdtUfCorr
(
    const VolField<scalar>& U,
    const SurfaceField<scalar>& Uf,
    const word
);

template<>
tmp<surfaceScalarField> boundedDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& U,
    const surfaceScalarField& phi,
    const word
);

template<>
tmp<surfaceScalarField> boundedDdtScheme<scalar>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& Uf,
    const word
);

template<>
tmp<surfaceScalarField> boundedDdtScheme<scalar>::fvcDdtPhiCorr
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
    #include "finiteVolume/ddtSchemes/boundedDdtScheme/boundedDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
