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
    Foam::wavePressureFvPatchScalarField

Description
    Boundary condition for the pressure gradient in a multiphase solver.
    To be used for wave boundary conditions. For more details see:

Description
    Boundary condition for the pressure in a multiphase solver.
    To be used for wave boundary conditions and is catered towards ship capsize
    simulations where the domain is rotating significantly. Boundary condition
    uses a user-specified threshold to switch between an inlet and atmosphere
    boundary condition as the mesh moves.

    The boundary condition is derived from mixedFvPatchField and it is generic
    in the sense, that it loads the abstract class waveTheory.

SourceFiles
    wavePressureCFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef wavePressureCFvPatchScalarField_H
#define wavePressureCFvPatchScalarField_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "fields/fvPatchFields/basic/mixed/mixedFvPatchFields.H"
#include "convexPolyhedral/convexPolyhedral.H"
#include "waveTheories/waveTheory/waveTheory.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class wavePressureCFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class wavePressureCFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public convexPolyhedral
{
    // Private data

        // Private objects
        autoPtr<waveTheories::waveTheory> waveProps_;

        //- Fixed Pressure
        scalarField p0_;

        //- Threshold for the dynamic boundary condition switch
        scalar thresholdValue;

        //- Threshold axis for the dynamic boundary condition switch (x,y,z)
        word thresholdAxis;


        // Private member functions

        //- Returns a scalarField of the signed distance to an arbitrary surface
        virtual void signedPointToSurfaceDistance
        (
            const pointField&,
            scalarField&
        );

        //- Returns a scalar of the signed distance to an arbitrary surface
        virtual scalar signedPointToSurfaceDistance
        (
            const point&
        ) const;

public:

    //- Runtime type information
    TypeName("wavePressureC");


    // Constructors

        //- Construct from patch and internal field
        wavePressureCFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wavePressureCFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given wavePressureCFvPatchScalarField onto
        //  a new patch
        wavePressureCFvPatchScalarField
        (
            const wavePressureCFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        wavePressureCFvPatchScalarField
        (
            const wavePressureCFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new wavePressureCFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        wavePressureCFvPatchScalarField
        (
            const wavePressureCFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new wavePressureCFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
