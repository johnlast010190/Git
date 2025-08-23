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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016-2024 Engys Ltd.

Class
    gaussFaLaplacianScheme

Description
    Basic second-order laplacian using face-gradients and Gauss' theorem.

SourceFiles
    gaussFaLaplacianScheme.C
    gaussFaLaplacianSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef gaussFaLaplacianScheme_H
#define gaussFaLaplacianScheme_H

#include "finiteArea/laplacianSchemes/faLaplacianScheme/faLaplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                       Class gaussLaplacianScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GType>
class gaussLaplacianScheme
:
    public fa::laplacianScheme<Type, GType>
{
    // Private Member Functions

        tmp<GeometricField<Type, faePatchField, faEdgeMesh>> gammaLnGradCorr
        (
                const edgeVectorField& SfGammaCorr,
                const GeometricField<Type, faPatchField, areaMesh>&
        );


public:

    //- Runtime type information
    TypeName("Gauss");


    // Constructors

        //- Construct null
        gaussLaplacianScheme(const faMesh& mesh)
        :
            laplacianScheme<Type, GType>(mesh)
        {}

        //- Construct from Istream
        gaussLaplacianScheme
        (
            const faMesh& mesh,
            Istream& is,
            const word& forceGradSchemeName
        )
        :
            laplacianScheme<Type, GType>(mesh, is, forceGradSchemeName)
        {}

        //- Construct from mesh, interpolation and snGradScheme schemes
        gaussLaplacianScheme
        (
            const faMesh& mesh,
            const tmp<edgeInterpolationScheme<GType>>& igs,
            const tmp<lnGradScheme<Type>>& lngs
        )
        :
            laplacianScheme<Type, GType>(mesh, igs, lngs)
        {}

        //- Disallow default bitwise copy construction
        gaussLaplacianScheme(const gaussLaplacianScheme&) = delete;


        //- Destructor
        virtual ~gaussLaplacianScheme()
        {}


    // Member Functions

        static tmp<faMatrix<Type>> famLaplacianUncorrected
        (
            const edgeScalarField& gammaMagSf,
            const edgeScalarField& deltaCoeffs,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<GeometricField<Type, faPatchField, areaMesh>> facLaplacian
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<faMatrix<Type>> famLaplacian
        (
            const GeometricField<GType, faePatchField, faEdgeMesh>&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<GeometricField<Type, faPatchField, areaMesh>> facLaplacian
        (
            const GeometricField<GType, faePatchField, faEdgeMesh>&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const gaussLaplacianScheme&) = delete;
};

// Use macros to emulate partial-specialisation of the the Laplacian functions
// for scalar diffusivity gamma

#define defineFamLaplacianScalarGamma(Type)                                 \
                                                                            \
template<>                                                                  \
tmp<faMatrix<Type>> gaussLaplacianScheme<Type, scalar>::famLaplacian       \
(                                                                           \
    const GeometricField<scalar, faePatchField, faEdgeMesh>&,              \
    const GeometricField<Type, faPatchField, areaMesh>&                      \
);                                                                          \
                                                                            \
template<>                                                                  \
tmp<GeometricField<Type, faPatchField, areaMesh>>                           \
gaussLaplacianScheme<Type, scalar>::facLaplacian                            \
(                                                                           \
    const GeometricField<scalar, faePatchField, faEdgeMesh>&,              \
    const GeometricField<Type, faPatchField, areaMesh>&                      \
);


defineFamLaplacianScalarGamma(scalar);
defineFamLaplacianScalarGamma(vector);
defineFamLaplacianScalarGamma(sphericalTensor);
defineFamLaplacianScalarGamma(symmTensor);
defineFamLaplacianScalarGamma(tensor);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "finiteArea/laplacianSchemes/gaussFaLaplacianScheme/gaussFaLaplacianScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
