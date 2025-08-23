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
    (c) 2024 Engys Ltd.

Class
    Foam::fv::LENCSnGrad

SourceFiles
    LENCSnGrad.C

Authors
    Simon Hill         <simon.hill@linde.com>
    Daniel Deising     <deising@mma.tu-darmstadt.de>
    All rights reserved.

Description
    Simple central-difference scalar-type snGrad scheme with implicit
    limited non-orthogonal and non-conjunctional correction.

\*---------------------------------------------------------------------------*/

#ifndef LENCSnGrad_H
#define LENCSnGrad_H

#include "finiteVolume/snGradSchemes/snGradScheme/snGradScheme.H"
// new
#include "interpolation/surfaceInterpolation/schemes/linear/linear.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "interpolation/surfaceInterpolation/schemes/skewCorrected/skewCorrectionVectors.H"
#include "dimensionedTypes/dimensionedType/dimensionedType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class LENCSnGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class LENCSnGrad
:
    public snGradScheme<Type>
{
    // Private Member Data

        //- Explicit Non-Orthogonal correction field
        // Has to be mutable that the const function deltaCoeffs
        // is allowed to change this variable
        mutable tmp<SurfaceField<Type>> tnonOrthCorr_;


public:

    //- Runtime type information
    TypeName("LENC");


    // Constructors

        //- Construct from mesh
        LENCSnGrad(const fvMesh& mesh)
        :
            snGradScheme<Type>(mesh),
            tnonOrthCorr_
            (
                surfaceScalarField::New
                (
                    "nonOrthCorr",
                    mesh,
                    dimensionSet(0, -1, 0, 0, 0, 0, 0)
                )
            )
        {}


        //- Construct from mesh and data stream
        LENCSnGrad
        (
            const fvMesh& mesh, Istream&, const word& forceGradSchemeName
        )
        :
            snGradScheme<Type>(mesh),
            tnonOrthCorr_
            (
                surfaceScalarField::New
                (
                    "nonOrthCorr",
                    mesh,
                    dimensionSet(0, -1, 0, 0, 0, 0, 0)
                )
            )
        {}


    //- Destructor
    virtual ~LENCSnGrad();

    // Member Functions

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> deltaCoeffs
        (
            const VolField<Type>& vf
        ) const;

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the LENCSnGrad
        //  for the given field using the gradients of the field components
        virtual tmp<SurfaceField<Type>>
        correction(const VolField<Type>&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LENCSnGrad&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "finiteVolume/snGradSchemes/LENCSnGrad/LENCSnGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
