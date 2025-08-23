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
    (c) 2024 Engys Ltd

Class
    Foam::fixedprghPressureFvPatchScalarField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    This boundary condition provides a boundary condition for the pressure (p)
    field that is equivalent to specifying 'fixedValue' on the p_rgh field.

Usage
    \table
        Property     | Description                | Required | Default value
        p_rgh        | Value of p_rgh             | yes      |
        hRef         | Reference height (at which p_rgh = p) | no | solver hRef
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fixedprghPressure;
        p_rgh           uniform 0;
        value           uniform 0; //Placeholder
    }
    \endverbatim


See also
    Foam::fixedValueFvPatchField

SourceFiles
    fixedprghPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedprghPressureFvPatchScalarField_H
#define fixedprghPressureFvPatchScalarField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class fixedprghPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fixedprghPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Specified value of p_rgh
        scalarField p_rgh_;

        //- hRef to override solver one, if specified in dictionary
        scalar hRef_;

        //- Whether hRef was specified in dictionary
        bool hRefSpecified_;


protected:

    // Protected functions

        Foam::tmp<Foam::scalarField> patchInternalRho();


public:

    //- Runtime type information
    TypeName("fixedprghPressure");


    // Constructors

        //- Construct from patch and internal field
        fixedprghPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedprghPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedprghPressureFvPatchScalarField
        //  onto a new patch
        fixedprghPressureFvPatchScalarField
        (
            const fixedprghPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedprghPressureFvPatchScalarField
        (
            const fixedprghPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const override
        {
            return tmp<fvPatchScalarField>
            (
                new fixedprghPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedprghPressureFvPatchScalarField
        (
            const fixedprghPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const override
        {
            return tmp<fvPatchScalarField>
            (
                new fixedprghPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return reference to the specified p_rgh pressure
            scalarField& p_rgh()
            {
                return p_rgh_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            ) override;

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            ) override;

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMapGIB(const gibFvPatchFieldMapper&) override;


        // Evaluation functions

            //- Inherit updateCoeffs from fixedValueFvPatchScalarField
            using fixedValueFvPatchScalarField::updateCoeffs;

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs() override;


        //- Write
        virtual void write(Ostream&) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
