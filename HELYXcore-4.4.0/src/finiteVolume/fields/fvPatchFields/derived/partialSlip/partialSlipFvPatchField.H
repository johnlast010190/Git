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
    (c) 2011-2017 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::partialSlipFvPatchField

Group
    grpWallBoundaryConditions grpGenericBoundaryConditions

Description
    This boundary condition provides a partial slip condition.  The amount of
    slip is controlled by a user-supplied field.

Usage
    \table
        Property      | Description             | Required    | Default value
        refValue      | Reference value at zero slip | no     | 0
        valueFraction | Fraction of refValue used for boundary [0-1] | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            partialSlip;
        refValue        uniform 0.001;
        valueFraction   uniform 0.1;
        value           uniform 0;
    }
    \endverbatim

See also
    Foam::transformFvPatchField

SourceFiles
    partialSlipFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef partialSlipFvPatchField_H
#define partialSlipFvPatchField_H

#include "fields/fvPatchFields/basic/transform/transformFvPatchField.H"
#include "fvMatrices/fvBlockMatrix/fvBlockMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class partialSlipFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class partialSlipFvPatchField
:
    public transformFvPatchField<Type>
{
    // Private data

        //- Reference value
        Field<Type> refValue_;

        //- Fraction (0-1) of value used for boundary condition
        scalarField valueFraction_;

        // Relax
        const scalar relax_;

public:

    //- Runtime type information
    TypeName("partialSlip");


    // Constructors

        //- Construct from patch and internal field
        partialSlipFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        partialSlipFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given partialSlipFvPatchField onto a new patch
        partialSlipFvPatchField
        (
            const partialSlipFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        partialSlipFvPatchField
        (
            const partialSlipFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const override
        {
            return tmp<fvPatchField<Type>>
            (
                new partialSlipFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        partialSlipFvPatchField
        (
            const partialSlipFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const override
        {
            return tmp<fvPatchField<Type>>
            (
                new partialSlipFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Attributes

            //- Return false: this patch field is not altered by assignment
            virtual bool assignable() const override
            {
                return false;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            ) override;

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            ) override;

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMapGIB(const gibFvPatchFieldMapper&) override;

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchField<Type>&) override;


        // Return defining fields

            virtual Field<Type>& refValue()
            {
                return refValue_;
            }

            virtual const Field<Type>& refValue() const
            {
                return refValue_;
            }

            virtual scalarField& valueFraction()
            {
                return valueFraction_;
            }

            virtual const scalarField& valueFraction() const
            {
                return valueFraction_;
            }


        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<Type>> snGrad() const override;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            ) override;

            //- Return face-gradient transform diagonal
            virtual tmp<Field<Type>> snGradTransformDiag() const override;

            //- Return face-gradient transform diagonal
            virtual tmp<Field<tensor>> transformationField() const override;

            //- Return the matrix diagonal coefficients
            virtual tmp<Field<Type>> valueDivInternalCoeffs
            (
                const tmp<scalarField>&
            ) const override;

            //- Return the matrix source coefficients
            virtual tmp<Field<Type>> valueDivBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const override;

        //- Relax the boundary condition
        virtual void boundaryRelaxMatrix(fvBlockMatrix<Type>&) const override {}


        //- Write
        virtual void write(Ostream&) const override;


    // Member operators

        virtual void operator=(const UList<Type>&) override {}

        virtual void operator=(const fvPatchField<Type>&) override {}
        virtual void operator+=(const fvPatchField<Type>&) override {}
        virtual void operator-=(const fvPatchField<Type>&) override {}
        virtual void operator*=(const fvPatchField<scalar>&) override {}
        virtual void operator/=(const fvPatchField<scalar>&) override {}

        virtual void operator+=(const Field<Type>&) override {}
        virtual void operator-=(const Field<Type>&) override {}

        virtual void operator*=(const Field<scalar>&) override {}
        virtual void operator/=(const Field<scalar>&) override {}

        virtual void operator=(const Type&) override {}
        virtual void operator+=(const Type&) override {}
        virtual void operator-=(const Type&) override {}
        virtual void operator*=(const scalar) override {}
        virtual void operator/=(const scalar) override {}
};

// * * * * * * * * * * * Template Specialisations  * * * * * * * * * * * * * //

template<>
void partialSlipFvPatchField<vector>::boundaryRelaxMatrix
(
    fvBlockMatrix<vector>&
) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
