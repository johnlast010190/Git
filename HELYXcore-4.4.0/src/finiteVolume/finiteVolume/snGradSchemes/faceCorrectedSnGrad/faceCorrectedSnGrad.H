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
    Foam::fv::faceCorrectedSnGrad

Group
    grpFvSnGradSchemes

Description
    Simple central-difference snGrad scheme with non-orthogonal correction.

SourceFiles
    faceCorrectedSnGrad.C

\*---------------------------------------------------------------------------*/

#ifndef faceCorrectedSnGrad_H
#define faceCorrectedSnGrad_H

#include "finiteVolume/snGradSchemes/snGradScheme/snGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class faceCorrectedSnGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class faceCorrectedSnGrad
:
    public snGradScheme<Type>
{
public:

    //- Runtime type information
    TypeName("faceCorrected");


    // Constructors

        //- Construct from mesh
        faceCorrectedSnGrad(const fvMesh& mesh)
        :
            snGradScheme<Type>(mesh)
        {}


        //- Construct from mesh and data stream
        faceCorrectedSnGrad
        (
            const fvMesh& mesh,
            Istream&,
            const word& forceGradSchemeName
        )
        :
            snGradScheme<Type>(mesh)
        {}


    //- Destructor
    virtual ~faceCorrectedSnGrad();


    // Member Functions

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> deltaCoeffs
        (
            const VolField<Type>&
        ) const
        {
            return this->mesh().nonOrthDeltaCoeffs();
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the faceCorrectedSnGrad
        //  for the given field using the gradient of the field
        tmp<SurfaceField<Type>>
        fullGradCorrection
        (
            const VolField<Type>&
        ) const;

        //- Return the explicit correction to the faceCorrectedSnGrad
        //  for the given field using the gradients of the field components
        virtual tmp<SurfaceField<Type>>
        correction(const VolField<Type>&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const faceCorrectedSnGrad&) = delete;
};


// * * * * * * * * Template Member Function Specialisations  * * * * * * * * //

template<>
tmp<surfaceScalarField> faceCorrectedSnGrad<scalar>::correction
(
    const volScalarField& vsf
) const;


template<>
tmp<surfaceVectorField> faceCorrectedSnGrad<vector>::correction
(
    const volVectorField& vvf
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/snGradSchemes/faceCorrectedSnGrad/faceCorrectedSnGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
