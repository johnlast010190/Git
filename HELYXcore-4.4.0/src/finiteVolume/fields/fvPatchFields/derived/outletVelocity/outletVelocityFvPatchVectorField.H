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
    Foam::outletVelocityFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    This boundary condition provides a velocity boundary condition at an outlet

Usage

    Example of the boundary condition specification
    \verbatim
    <patchName>
    {
        type                outletVelocity;
        value               uniform (0 0 0);
        relax               0.2;
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField

SourceFiles
    outletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef outletVelocityFvPatchVectorField_H
#define outletVelocityFvPatchVectorField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "fields/fvPatchFields/derived/flowRateBase/flowRateBase.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class outletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class outletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- patch relax in block solver
        scalar relax_;

public:

   //- Runtime type information
   TypeName("outletVelocity");


   // Constructors

        //- Construct from patch and internal field
        outletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        outletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  outletVelocityFvPatchVectorField
        //  onto a new patch
        outletVelocityFvPatchVectorField
        (
            const outletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        outletVelocityFvPatchVectorField
        (
            const outletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const override
        {
            return tmp<fvPatchVectorField>
            (
                new outletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        outletVelocityFvPatchVectorField
        (
            const outletVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const override
        {
            return tmp<fvPatchVectorField>
            (
                new outletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap
        (
            const fvPatchFieldMapper&
        ) override;

        //- Reverse map the given fvPatchField onto this fvPatchField
        virtual void rmap
        (
            const fvPatchVectorField&,
            const labelList&
        ) override;

        //- Reset the fvPatchField to the given fvPatchField
        //  Used for mesh to mesh mapping
        virtual void reset(const fvPatchVectorField&) override;

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMapGIB(const gibFvPatchFieldMapper&) override;


        //- Return the matrix diagonal coefficients
        virtual tmp<Field<vector>> valueInternalCoeffs
        (
            const tmp<scalarField>&
        ) const override;

        //- Return the matrix source coefficients
        virtual tmp<Field<vector>> valueBoundaryCoeffs
        (
            const tmp<scalarField>&
        ) const override;

        //- Return the matrix diagonal coefficients
        virtual tmp<Field<vector>> valueDivInternalCoeffs
        (
            const tmp<scalarField>&
        ) const override;

        //- Return the matrix source coefficients
        virtual tmp<Field<vector>> valueDivBoundaryCoeffs
        (
            const tmp<scalarField>&
        ) const override;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs() override;

        //- Manipulate matrix with given weights
        virtual void boundaryRelaxMatrix(fvBlockMatrix<vector>&) const override;

        //- Write
        virtual void write(Ostream&) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
