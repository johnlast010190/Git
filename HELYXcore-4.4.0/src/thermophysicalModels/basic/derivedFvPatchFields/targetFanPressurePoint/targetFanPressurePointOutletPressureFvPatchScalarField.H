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
    (c) 2020-2024 Engys Ltd.

Class
    Foam::targetFanPressurePointOutletPressureFvPatchScalarField

Description
    Foam::targetFanPressurePointOutletPressureFvPatchScalarField

    Modify pressure to satisfy user-specified pressure
    at a reference point.

Usage

    Example of the boundary condition specification for a ref. point:
    \verbatim
    myPatch
    {
        type            targetFanPressurePoint;
        initialPressure -28000;
        referencePoint  ( 0.5 1.8 0 );
        timeFanCurve    false;
        fanCurve        table
            (
                ( 0 400 )
                ( 2.5 420 )
                ( 5 380 )
                ( 7.5 300 )
                ( 9 0 )
            );
        outOfBounds     clamp;
    }
    \endverbatim


SourceFiles
    targetFanPressurePointOutletPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef targetFanPressurePointOutletPressureFvPatchScalarField_H
#define targetFanPressurePointOutletPressureFvPatchScalarField_H

#include "fields/fvPatchFields/basic/mixed/mixedFvPatchFields.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "primitives/functions/Function2/Function2/Function2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class targetFanPressurePointOutletPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class targetFanPressurePointOutletPressureFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Name of the flux
        word phiName_;

        word rhoName_;

        const point refPoint_;

        scalar initialP_;

        //- Is the fan curve dependent on time
        Switch isTimeCurve_;

        //- Tabulated fan curve
        autoPtr<Function1<scalar>> fanCurve_;

        //- Tabulated dependent on time/flow rate
        autoPtr<Function2<scalar>> timeFanCurve_;

        bool isInitial_;

    // Private member functions

        //- Function checks flux dimensions and returns if
        //  they are massflux
        bool checkPhiDimensions
        (
            const fvsPatchField<scalar>& phiP
        ) const;

        scalar referenceValue
        (
            const word& fieldName
        );

//        static const label closestCell(const word& fieldName);

public:

    //- Runtime type information
    TypeName("targetFanPressurePoint");


    // Constructors

        //- Construct from patch and internal field
        targetFanPressurePointOutletPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        targetFanPressurePointOutletPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  targetFanPressurePointOutletPressureFvPatchScalarField onto a new patch
        targetFanPressurePointOutletPressureFvPatchScalarField
        (
            const targetFanPressurePointOutletPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        targetFanPressurePointOutletPressureFvPatchScalarField
        (
            const targetFanPressurePointOutletPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new targetFanPressurePointOutletPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        targetFanPressurePointOutletPressureFvPatchScalarField
        (
            const targetFanPressurePointOutletPressureFvPatchScalarField&,
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
                new targetFanPressurePointOutletPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

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
