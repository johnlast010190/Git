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
    (c) 2022-2024 Engys Ltd

Class
    Foam::multivariateDeferred

Group
    grpFvmultivariateDeferreds

Description
    Deferred-correction scheme - turns any interpolation scheme into an upwind-
    weighted scheme with and explicit multivariateDeferred correction.

SourceFiles
    multivariateDeferred.C

\*---------------------------------------------------------------------------*/

#ifndef multivariateDeferred_H
#define multivariateDeferred_H

#include "interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.H"
#include "interpolation/surfaceInterpolation/multivariateSchemes/upwind/multivariateUpwind.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class multivariateDeferred Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class multivariateDeferred
:
    public multivariateUpwind<Type>
{
protected:

    // Protected data

        tmp<multivariateSurfaceInterpolationScheme<Type>> scheme_;


public:

    //- Runtime type information
    TypeName("deferred");


    // Constructors

        //- Construct for field, faceFlux and Istream
        multivariateDeferred
        (
            const fvMesh& mesh,
            const typename multivariateSurfaceInterpolationScheme<Type>::
                fieldTable& fields,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        )
        :
            multivariateUpwind<Type>
            (
                mesh,
                fields,
                faceFlux,
                schemeData
            ),
            scheme_
            (
                multivariateSurfaceInterpolationScheme<Type>::New
                (
                    mesh, fields, faceFlux, schemeData
                )
            )
        {}

        //- Disallow default bitwise copy construction
        multivariateDeferred(const multivariateDeferred&) = delete;


    //- Destructor
    virtual ~multivariateDeferred() {};


    // Member operators

        //- surfaceInterpolationScheme sub-class returned by operator(field)
        //  which is used as the interpolation scheme for the field
        class fieldScheme
        :
            public multivariateUpwind<Type>::fieldScheme
        {
            // Private data

                tmp<surfaceInterpolationScheme<Type>> scheme_;

        public:

            // Constructors

                fieldScheme
                (
                    const VolField<Type>& field,
                    const surfaceScalarField& faceFlux,
                    const multivariateSurfaceInterpolationScheme<Type>& scheme
                )
                :
                    multivariateUpwind<Type>::fieldScheme(field, faceFlux),
                    scheme_(scheme(field))
                {}


            // Member Functions

                //- Return true if this scheme uses an explicit correction
                virtual bool corrected() const override
                {
                    return true;
                }

                //- Return explicit correction
                tmp<SurfaceField<Type>> correction
                (
                    const VolField<Type>& vf
                ) const override;
        };

        tmp<surfaceInterpolationScheme<Type>> operator()
        (
            const VolField<Type>& field
        ) const
        {
            return tmp<surfaceInterpolationScheme<Type>>
            (
                new fieldScheme(field, this->faceFlux_, scheme_())
            );
        }

        //- Disallow default bitwise assignment
        void operator=(const multivariateDeferred&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interpolation/surfaceInterpolation/multivariateSchemes/deferred/multivariateDeferredTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
