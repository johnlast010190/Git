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
    (c) 2011 OpenFOAM Foundation
    (c) 2010-2024 Engys Ltd.

Class
    Foam::omegaGaussLaplacianScheme

Description
    omegaGaussLaplacianScheme interpolation scheme class derived from upwind and returns
    upwind weighting factors and also applies a gradient-based explicit
    correction.

SourceFiles
    omegaGaussLaplacianScheme.C

\*---------------------------------------------------------------------------*/

#ifndef omegaGaussLaplacianScheme_H
#define omegaGaussLaplacianScheme_H

#include "finiteVolume/laplacianSchemes/laplacianScheme/laplacianScheme.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{
/*---------------------------------------------------------------------------*\
                Class omegaGaussLaplacianScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GType>
class omegaGaussLaplacianScheme
:
public fv::laplacianScheme<Type, GType>
{

    // Private Member Functions

        //- Calculate the gammaSnGradCorr field
        tmp<SurfaceField<Type>> gammaSnGradCorr
        (
            const surfaceVectorField& SfGammaCorr,
            const VolField<Type>&
        );

        scalar getDeltaStarCoefficient
        (
            const label& patchID,
            const label& patchFaceLabel,
            const label& faceLabel,
            const scalar& yPlus
        );


public:

    //- Runtime type information
    TypeName("omegaGauss");


    // Constructors

        //- Construct null
        omegaGaussLaplacianScheme(const fvMesh& mesh)
        :
            laplacianScheme<Type, GType>(mesh)
        {}

        //- Construct from faceFlux
        omegaGaussLaplacianScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux
        )
        :
            laplacianScheme<Type, GType>(mesh)
        {}

        //- Construct from Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        omegaGaussLaplacianScheme
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& schemeData,
            const word& forceGradSchemeName=word::null
        )
        :
            laplacianScheme<Type, GType>(mesh)
        {}

        //- Construct from mesh, interpolation and snGradScheme schemes
        omegaGaussLaplacianScheme
        (
            const fvMesh& mesh,
            const tmp<surfaceInterpolationScheme<GType>>& igs,
            const tmp<snGradScheme<Type>>& sngs
        )
        :
            laplacianScheme<Type, GType>(mesh, igs, sngs)
        {}

        //- Disallow default bitwise copy construction
        omegaGaussLaplacianScheme(const omegaGaussLaplacianScheme&) = delete;


    //- Destructor
    virtual ~omegaGaussLaplacianScheme()
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

        tmp<fvMatrix<Type>> fvmLaplacian
        (
            const SurfaceField<GType>&,
            const volScalarField&,
            const VolField<Type>&
        )
        {
            NotImplemented;
        }

        tmp<fvBlockMatrix<Type>> fvmBLaplacian
        (
            const SurfaceField<GType>&,
            const VolField<Type>&
        )
        {
            FatalErrorInFunction
                << "Not Implemented " << abort(FatalError);
            return tmp<fvBlockMatrix<Type>>();
        }

        tmp<VolField<Type>> fvcLaplacian
        (
            const SurfaceField<GType>&,
            const VolField<Type>&
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const omegaGaussLaplacianScheme&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "omegaGaussLaplacianScheme.C"
#endif

#endif

// ************************************************************************* //
