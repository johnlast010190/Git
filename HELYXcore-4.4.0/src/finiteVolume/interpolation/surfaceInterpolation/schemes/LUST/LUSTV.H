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
    (c) 2017-2024 Engys Ltd.

Class
    Foam::LUSTV

Description
    LUSTV: Linear-upwindV stabilised transport.

    Interpolation scheme class derived from linearUpwindV which returns blended
    linear/linear-upwindV weighting factors and also applies a explicit
    gradient-based correction obtained from the linearUpwind scheme.  The
    blending-factor is set to 0.75 linear which optimises the balance between
    accuracy and stability on a range of LES cases with a range of mesh quality.

SourceFiles
    LUSTV.C

\*---------------------------------------------------------------------------*/

#ifndef LUSTV_H
#define LUSTV_H

#include "interpolation/surfaceInterpolation/schemes/linearUpwind/linearUpwindV.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class LUST Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class LUSTV
:
    public linearUpwindV<Type>
{
public:

    //- Runtime type information
    TypeName("LUSTV");


    // Constructors

        //- Construct from mesh and Istream
        LUSTV
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& schemeData
        )
        :
            linearUpwindV<Type>(mesh, db, schemeData)
        {}

        //- Construct from mesh, faceFlux and Istream
        LUSTV
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        )
        :
            linearUpwindV<Type>(mesh, faceFlux, schemeData)
        {}

        //- Disallow default bitwise copy construction
        LUSTV(const LUSTV&) = delete;


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> weights
        (
            const VolField<Type>&
        ) const
        {
            return
                0.75*this->mesh().surfaceInterpolation::weights()
              + 0.25*linearUpwindV<Type>::weights();
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<SurfaceField<Type>>
        correction
        (
            const VolField<Type>& vf
        ) const
        {
            return 0.25*linearUpwindV<Type>::correction(vf);
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LUSTV&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
