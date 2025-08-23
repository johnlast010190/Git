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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2020 Engys Ltd.

Class
    Foam::basicSymmetryFvPatchField

Group
    grpConstraintBoundaryConditions

Description
    A symmetry patch

SourceFiles
    basicSymmetryFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef basicSymmetryFvPatchField_H
#define basicSymmetryFvPatchField_H

#include "fields/fvPatchFields/basic/transform/transformFvPatchField.H"
#include "fvMesh/fvPatches/constraint/symmetry/symmetryFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class basicSymmetryFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class basicSymmetryFvPatchField
:
    public transformFvPatchField<Type>
{

public:

    // Constructors

        //- Construct from patch and internal field
        basicSymmetryFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        basicSymmetryFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given basicSymmetryFvPatchField onto a new
        //  patch
        basicSymmetryFvPatchField
        (
            const basicSymmetryFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        basicSymmetryFvPatchField
        (
            const basicSymmetryFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const override
        {
            return tmp<fvPatchField<Type>>
            (
                new basicSymmetryFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        basicSymmetryFvPatchField
        (
            const basicSymmetryFvPatchField<Type>&,
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
                new basicSymmetryFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        //- Return gradient at boundary
        virtual tmp<Field<Type>> snGrad() const override;

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
        ) override;

        //- Return face-gradient transform diagonal
        virtual tmp<Field<Type>> snGradTransformDiag() const override;

        //- Return face-gradient transform diagonal
        virtual tmp<Field<tensor>> transformationField() const override;

        //- Return patch value for explicit gradient evaluation
        tmp<Field<Type>> gradientBoundaryValue() const override;
};


// * * * * * * * * * * * Template Specialisations  * * * * * * * * * * * * * //

template<>
tmp<scalarField > basicSymmetryFvPatchField<scalar>::snGrad() const;

template<>
void basicSymmetryFvPatchField<scalar>::evaluate
(
    const Pstream::commsTypes commsType
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fields/fvPatchFields/basic/basicSymmetry/basicSymmetryFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
