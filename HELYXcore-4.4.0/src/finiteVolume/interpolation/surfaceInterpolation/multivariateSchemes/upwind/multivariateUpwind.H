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
    (c) 2022-2024 Engys Ltd.

Class
    Foam::multivariateUpwind

Description
    Upwind-difference form of the multivariate surfaceInterpolationScheme.

SourceFiles
    multivariateUpwindmake.C

\*---------------------------------------------------------------------------*/

#ifndef multivariateUpwind_H
#define multivariateUpwind_H

#include "interpolation/surfaceInterpolation/multivariateSchemes/multivariateSurfaceInterpolationScheme/multivariateSurfaceInterpolationScheme.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class multivariateUpwind Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class multivariateUpwind
:
    public multivariateSurfaceInterpolationScheme<Type>
{
protected:

    // Protected data

        const surfaceScalarField& faceFlux_;


public:

    //- Runtime type information
    TypeName("upwind");


    // Constructors

        //- Construct for field, faceFlux and Istream
        multivariateUpwind
        (
            const fvMesh& mesh,
            const typename multivariateSurfaceInterpolationScheme<Type>::
                fieldTable& fields,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        )
        :
            multivariateSurfaceInterpolationScheme<Type>
            (
                mesh,
                fields,
                faceFlux,
                schemeData
            ),
            faceFlux_(faceFlux)
        {}

        //- Disallow default bitwise copy construction
        multivariateUpwind(const multivariateUpwind&) = delete;


    // Member Operators

        //- surfaceInterpolationScheme sub-class returned by operator(field)
        //  which is used as the interpolation scheme for the field
        class fieldScheme
        :
            public multivariateSurfaceInterpolationScheme<Type>::fieldScheme
        {
        protected:

            // Protected data

                const surfaceScalarField& faceFlux_;

        public:

            // Constructors

                fieldScheme
                (
                    const VolField<Type>& field,
                    const surfaceScalarField& faceFlux
                )
                :
                    multivariateSurfaceInterpolationScheme<Type>::
                        fieldScheme(field),
                    faceFlux_(faceFlux)
                {}


            // Member Functions

                //- Return the interpolation weighting factors
                tmp<surfaceScalarField> weights
                (
                    const VolField<Type>&
                ) const
                {
                    return pos0(faceFlux_);
                }
        };

        tmp<surfaceInterpolationScheme<Type>> operator()
        (
            const VolField<Type>& field
        ) const
        {
            return tmp<surfaceInterpolationScheme<Type>>
            (
                new fieldScheme(field, faceFlux_)
            );
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multivariateUpwind&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
