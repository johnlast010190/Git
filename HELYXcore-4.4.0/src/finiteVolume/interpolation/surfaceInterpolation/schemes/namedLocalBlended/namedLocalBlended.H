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
    (c) 2013-2024 Engys Ltd.

Class
    Foam::namedLocalBlended

Description
    Two-scheme namedLocalBlended differencing scheme.

SourceFiles
    namedLocalBlended.C

\*---------------------------------------------------------------------------*/

#ifndef namedLocalBlended_H
#define namedLocalBlended_H

#include "interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class namedLocalBlended Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class namedLocalBlended
:
    public surfaceInterpolationScheme<Type>
{
    // Private Member Functions

        word blendingFieldName_;

        //- Scheme 1
        tmp<surfaceInterpolationScheme<Type>> tScheme1_;

        //- Scheme 2
        tmp<surfaceInterpolationScheme<Type>> tScheme2_;


public:

    //- Runtime type information
    TypeName("namedLocalBlended");


    // Constructors

        //- Construct from mesh and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        namedLocalBlended
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            blendingFieldName_(is),
            tScheme1_
            (
                surfaceInterpolationScheme<Type>::New(mesh, db, is)
            ),
            tScheme2_
            (
                surfaceInterpolationScheme<Type>::New(mesh, db, is)
            )
        {}

        //- Construct from mesh, faceFlux and Istream
        namedLocalBlended
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            blendingFieldName_(is),
            tScheme1_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            ),
            tScheme2_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            )
        {}

        //- Disallow default bitwise copy construction
        namedLocalBlended(const namedLocalBlended&) = delete;


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const VolField<Type>& vf
        ) const
        {
            const surfaceScalarField& blendingFactor =
                this->mesh().objectRegistry::template
                lookupObject<const surfaceScalarField>
                (
                    blendingFieldName_
                );

            return
                blendingFactor*tScheme1_().weights(vf)
              + (scalar(1) - blendingFactor)*tScheme2_().weights(vf);
        }

        //- Return the face-interpolate of the given cell field
        //  with explicit correction
        tmp<SurfaceField<Type>>
        interpolate(const VolField<Type>& vf) const
        {
            const surfaceScalarField& blendingFactor =
            (
                this->mesh().objectRegistry::template
                lookupObject<const surfaceScalarField>
                (
                    blendingFieldName_
                )
            );

            return
                blendingFactor*tScheme1_().interpolate(vf)
              + (scalar(1) - blendingFactor)*tScheme2_().interpolate(vf);
        }


        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return tScheme1_().corrected() || tScheme2_().corrected();
        }


        //- Return the explicit correction to the face-interpolate
        //  for the given field
        virtual tmp<SurfaceField<Type>>
        correction
        (
            const VolField<Type>& vf
        ) const
        {
            const surfaceScalarField& blendingFactor =
                this->mesh().objectRegistry::template
                lookupObject<const surfaceScalarField>
                (
                    blendingFieldName_
                );

            if (tScheme1_().corrected())
            {
                if (tScheme2_().corrected())
                {
                    return
                    (
                        blendingFactor
                      * tScheme1_().correction(vf)
                      + (scalar(1.0) - blendingFactor)
                      * tScheme2_().correction(vf)
                    );
                }
                else
                {
                    return
                    (
                        blendingFactor
                      * tScheme1_().correction(vf)
                    );
                }
            }
            else if (tScheme2_().corrected())
            {
                return
                (
                    (scalar(1.0) - blendingFactor)
                  * tScheme2_().correction(vf)
                );
            }
            else
            {
                return tmp<SurfaceField<Type>>
                (
                    nullptr
                );
            }
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const namedLocalBlended&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
