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
    (c) 2017-2024 Engys Ltd.

Class
    Foam::fv::phasicGaussGrad

Group
    grpFvGradSchemes

Description
    Basic second-order gradient scheme using face-interpolation
    and Gauss' theorem.

SourceFiles
    phasicGaussGrad.C

\*---------------------------------------------------------------------------*/

#ifndef phasicGaussGrad_H
#define phasicGaussGrad_H

#include "finiteVolume/gradSchemes/gradScheme/gradScheme.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.H"
#include "interpolation/surfaceInterpolation/schemes/linear/linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class phasicGaussGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class phasicGaussGrad
:
    public fv::gradScheme<Type>
{
    // Private data

        tmp<surfaceInterpolationScheme<Type>> tinterpScheme_;

        //minimum value of phase fraction for which field operations are
        //defined
        scalar alphaCutoff_;

        //phase fraction field name
        word alphaName_;


    // Private Member Functions

        //- Calculate and set phase space weights
        static void phasicFaceWeights
        (
            const scalar,
            const scalar,
            const scalar,
            scalar&,
            scalar&,
            scalar&
        ) ;


public:

    //- Runtime type information
    TypeName("phasicGauss");


    // Constructors

        //- Construct from mesh
        phasicGaussGrad(const fvMesh& mesh)
        :
            gradScheme<Type>(mesh),
            tinterpScheme_(new linear<Type>(mesh)),
            alphaCutoff_(1e-5),
            alphaName_("alpha.foam")
        {
            NotImplemented;
        }

        //- Construct from mesh and Istream
        phasicGaussGrad
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& is
        )
        :
            gradScheme<Type>(mesh),
            tinterpScheme_(nullptr),
            alphaCutoff_(readScalar(is)),
            alphaName_(is)
        {
            if (is.eof())
            {
                tinterpScheme_ =
                    tmp<surfaceInterpolationScheme<Type>>
                    (
                        new linear<Type>(mesh)
                    );
            }
            else
            {
                tinterpScheme_ =
                    tmp<surfaceInterpolationScheme<Type>>
                    (
                        surfaceInterpolationScheme<Type>::New(mesh, db, is)
                    );
            }
        }

        //- Disallow default bitwise copy construction
        phasicGaussGrad(const phasicGaussGrad&) = delete;


    // Member Functions

        //- Surface field internal face gradient contribution
        static tmp<VolField<typename outerProduct<vector, Type>::type>>
        gradFinternal
        (
            const SurfaceField<Type>&,
            const VolField<Type>&,
            const word& name,
            const scalar alphaCutoff,
            const word alphaName
        );

        //- Return the gradient of the given field
        //  calculated using Gauss' theorem on the given surface field
        static tmp<VolField<typename outerProduct<vector, Type>::type>> gradf
        (
            const SurfaceField<Type>&,
            const VolField<Type>&,
            const word& name
        );

        //- Return the gradient of the given field to the gradScheme::grad
        //  for optional caching
        virtual tmp<VolField<typename outerProduct<vector, Type>::type>>
        calcGrad
        (
            const VolField<Type>& vsf,
            const word& name
        ) const;

        //- Correct the boundary values of the gradient using the patchField
        // snGrad functions
        static void correctBoundaryConditions
        (
            const VolField<Type>&,
            VolField<typename outerProduct<vector, Type>::type>&
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const phasicGaussGrad&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/gradSchemes/phasicGaussGrad/phasicGaussGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
