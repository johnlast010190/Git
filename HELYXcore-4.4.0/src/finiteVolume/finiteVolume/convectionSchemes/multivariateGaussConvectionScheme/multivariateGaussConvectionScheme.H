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

Class
    Foam::fv::multivariateGaussConvectionScheme

Group
    grpFvConvectionSchemes

Description
    Basic second-order convection using face-gradients and Gauss' theorem.

SourceFiles
    multivariateGaussConvectionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef multivariateGaussConvectionScheme_H
#define multivariateGaussConvectionScheme_H

#include "finiteVolume/convectionSchemes/convectionScheme/convectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class multivariateGaussConvectionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class multivariateGaussConvectionScheme
:
    public fv::convectionScheme<Type>
{
    // Private data

        tmp<multivariateSurfaceInterpolationScheme<Type>> tinterpScheme_;


public:

    //- Runtime type information
    TypeName("Gauss");


    // Constructors

        //- Construct from flux and Istream
        multivariateGaussConvectionScheme
        (
            const fvMesh& mesh,
            const typename multivariateSurfaceInterpolationScheme<Type>::
                fieldTable& fields,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            tinterpScheme_
            (
                multivariateSurfaceInterpolationScheme<Type>::New
                (
                    mesh, fields, faceFlux, is
                )
            )
        {}


    // Member Functions

        tmp<multivariateSurfaceInterpolationScheme<Type>>
        interpolationScheme() const
        {
            return tinterpScheme_;
        }

        tmp<SurfaceField<Type>> interpolate
        (
            const surfaceScalarField&,
            const VolField<Type>&
        ) const;

        tmp<SurfaceField<Type>> flux
        (
            const surfaceScalarField&,
            const VolField<Type>&
        ) const;

        tmp<fvMatrix<Type>> fvmDiv
        (
            const surfaceScalarField&,
            const VolField<Type>&
        ) const;

        tmp<fvMatrix<Type>> fvmDiv
        (
            const surfaceScalarField&,
            const volScalarField&,
            const VolField<Type>&
        ) const;

        tmp<VolField<Type>> fvcDiv
        (
            const surfaceScalarField&,
            const VolField<Type>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
