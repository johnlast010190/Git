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
    (c) 2019-2024 Engys Ltd.

Class
    Foam::fv::compactGaussLaplacianScheme

Group
    grpFvLaplacianSchemes

Description
    (An)isotropic Laplacian using compact stencil formulation with over-
    relaxation for face non-orthogonality (if 'corrected' snGrad scheme is
    used).

SourceFiles
    compactGaussLaplacianScheme.C

\*---------------------------------------------------------------------------*/

#ifndef compactGaussLaplacianScheme_H
#define compactGaussLaplacianScheme_H

#include "finiteVolume/laplacianSchemes/laplacianScheme/laplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                Class compactGaussLaplacianScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GType>
class compactGaussLaplacianScheme
:
    public fv::laplacianScheme<Type, GType>
{
public:

    //- Runtime type information
    TypeName("compactGauss");


    // Constructors

        //- Construct null
        compactGaussLaplacianScheme(const fvMesh& mesh)
        :
            laplacianScheme<Type, GType>(mesh)
        {}

        //- Construct from Istream
        compactGaussLaplacianScheme
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& is,
            const word& forceGradSchemeName=word::null
        )
        :
            laplacianScheme<Type, GType>(mesh, db, is, forceGradSchemeName)
        {}

        //- Construct from mesh, interpolation and snGradScheme schemes
        compactGaussLaplacianScheme
        (
            const fvMesh& mesh,
            const tmp<surfaceInterpolationScheme<GType>>& igs,
            const tmp<snGradScheme<Type>>& sngs
        )
        :
            laplacianScheme<Type, GType>(mesh, igs, sngs)
        {}

        //- Disallow default bitwise copy construction
        compactGaussLaplacianScheme
        (
            const compactGaussLaplacianScheme&
        ) = delete;


    //- Destructor
    virtual ~compactGaussLaplacianScheme()
    {}


    // Member Functions

        tmp<VolField<Type>> fvcLaplacian
        (
            const VolField<Type>&
        );

        tmp<fvMatrix<Type>> fvmLaplacian
        (
            const SurfaceField<GType>&,
            const VolField<Type>&
        );

        tmp
        <
            fvBlockMatrix<Type>
        > fvmBLaplacian
        (
            const SurfaceField<GType>&,
            const VolField<Type>&
        )
        {
            FatalErrorInFunction << "Not Implemented " << abort(FatalError);
            return tmp<fvBlockMatrix<Type>>();
        }

        tmp<VolField<Type>> fvcLaplacian
        (
            const SurfaceField<GType>&,
            const VolField<Type>&
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const compactGaussLaplacianScheme&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "compactGaussLaplacianScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
