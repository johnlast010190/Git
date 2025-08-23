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
    (c) 2024 Engys Ltd.

Class
    Foam::linearUpwind

Group
    grpFvSurfaceInterpolationSchemes

Description
    linearUpwind interpolation scheme class derived from upwind and returns
    upwind weighting factors and also applies a gradient-based explicit
    correction.

SourceFiles
    linearUpwind.C

\*---------------------------------------------------------------------------*/

#ifndef linearUpwind_H
#define linearUpwind_H

#include "interpolation/surfaceInterpolation/limitedSchemes/upwind/upwind.H"
#include "finiteVolume/gradSchemes/gaussGrad/gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class linearUpwind Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class linearUpwind
:
    public upwind<Type>
{
    // Private Data

        word gradSchemeName_;


public:

    //- Runtime type information
    TypeName("linearUpwind");


    // Constructors

        //- Construct from faceFlux
        linearUpwind
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux
        )
        :
            upwind<Type>(mesh, faceFlux),
            gradSchemeName_("grad")
        {}

        //- Construct from Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        linearUpwind
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& schemeData
        )
        :
            upwind<Type>(mesh, db, schemeData),
            gradSchemeName_(schemeData)
        {}

        //- Construct from faceFlux and Istream
        linearUpwind
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        )
        :
            upwind<Type>(mesh, faceFlux, schemeData),
            gradSchemeName_(schemeData)
        {}

        //- Disallow default bitwise copy construction
        linearUpwind(const linearUpwind&) = delete;


    // Member Functions

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<SurfaceField<Type>>
        correction
        (
            const VolField<Type>&
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const linearUpwind&) = delete;
};


// Specialize for volVectorField to support cached gradients
template<>
tmp<surfaceVectorField> linearUpwind<vector>::correction
(
    const volVectorField& vf
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
