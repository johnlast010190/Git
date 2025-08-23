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
    (c) 2012-2016 OpenFOAM Foundation
    (c) 2022-2025 Engys Ltd.

Class
    Foam::fv::scaledDdtScheme

Group
    grpFvDdtSchemes

Description
    scaled temporal ddt scheme.

    This allows to some equations to be solved using a different timestep in
    temporal term. The scale allows to use a temporal scheme with timestep
    value of scale*dt, where dt is selected from controlDict.

Usage

    Example of the scaled temporal scheme:
    \verbatim
        ddt(rho,K) scaled Euler 100;
    \endverbatim


SourceFiles
    scaledDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef scaledDdtScheme_H
#define scaledDdtScheme_H

#include "finiteVolume/ddtSchemes/ddtScheme/ddtScheme.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class scaledDdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class scaledDdtScheme
:
    public fv::ddtScheme<Type>
{
    // Private data

        tmp<fv::ddtScheme<Type>> scheme_;

        autoPtr<Function1<scalar>> scale_;


public:

    //- Runtime type information
    TypeName("scaled");


    // Constructors

        //- Construct from mesh and Istream
        scaledDdtScheme(const fvMesh& mesh, Istream& is);

        //- Disallow default bitwise copy construction
        scaledDdtScheme(const scaledDdtScheme&) = delete;


    // Member Functions

        scalar scale() const
        {
            const scalar curScale = scale_->value(mesh().time().value());
            scalar rScalar = 0.0;
            if (curScale!=0)
            {
                rScalar = 1.0/curScale;
            }
            return rScalar;
        }

        virtual scalar pseudoDeltaTValue() const
        {
            const scalar deltaT = mesh().time().deltaTValue();
            const scalar curScale = scale_->value(mesh().time().value());
            return curScale*deltaT;
        }

        //- Return mesh reference
        virtual const fvMesh& mesh() const
        {
            return fv::ddtScheme<Type>::mesh();
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
        void operator=(const scaledDdtScheme&) = delete;
};


template<>
tmp<surfaceScalarField> scaledDdtScheme<scalar>::fvcDdtUfCorr
(
    const VolField<scalar>& U,
    const SurfaceField<scalar>& Uf,
    const word
);

template<>
tmp<surfaceScalarField> scaledDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& U,
    const surfaceScalarField& phi,
    const word
);

template<>
tmp<surfaceScalarField> scaledDdtScheme<scalar>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& Uf,
    const word
);

template<>
tmp<surfaceScalarField> scaledDdtScheme<scalar>::fvcDdtPhiCorr
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
    #include "finiteVolume/ddtSchemes/scaledDdtScheme/scaledDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
