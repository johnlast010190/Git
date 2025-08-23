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
    (c) 2024 Engys Ltd.

Class
    Foam::nonConformalDiscreteMixingFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition enforces a discrete-mixing condition between a pair
    of non-conformal boundaries.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            nonConformalDiscreteMixing;
    }
    \endverbatim

Note
    The patches are allowed to be topologically dissimilar, i.e. if the owner
    patch is transformed to the neighbour patch, the patches do not necessarily
    need to be identical.

SourceFiles
    nonConformalDiscreteMixingFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalDiscreteMixingFvPatchField_H
#define nonConformalDiscreteMixingFvPatchField_H

#include "fields/fvPatchFields/basic/coupled/coupledFvPatchField.H"
#include "fvMesh/fvPatches/constraint/nonConformalDiscreteMixing/nonConformalDiscreteMixingFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class nonConformalDiscreteMixingFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class nonConformalDiscreteMixingFvPatchField
:
    public coupledFvPatchField<Type>
{
    // Private data

        //- Local reference cast into the non-conformal discrete-mixing patch
        const nonConformalDiscreteMixingFvPatch&
            nonConformalDiscreteMixingPatch_;


public:

    //- Runtime type information
    TypeName(nonConformalDiscreteMixingFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        nonConformalDiscreteMixingFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary.
        nonConformalDiscreteMixingFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given nonConformalDiscreteMixingFvPatchField
        //  onto a new patch
        nonConformalDiscreteMixingFvPatchField
        (
            const nonConformalDiscreteMixingFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nonConformalDiscreteMixingFvPatchField
        (
            const nonConformalDiscreteMixingFvPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        nonConformalDiscreteMixingFvPatchField
        (
            const nonConformalDiscreteMixingFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const override
        {
            return tmp<fvPatchField<Type>>
            (
                new nonConformalDiscreteMixingFvPatchField<Type>(*this)
            );
        }

        //- Construct and return a clone, setting internal field reference.
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const override
        {
            return tmp<fvPatchField<Type>>
            (
                new nonConformalDiscreteMixingFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return local reference cast into the non-conformal patch
            const nonConformalDiscreteMixingFvPatch& nonConformalPatch() const
            {
                return nonConformalDiscreteMixingPatch_;
            }


        // Evaluation functions

            //- Return neighbour coupled internal cell data
            virtual tmp<Field<Type>> patchNeighbourField
            (
                const Pstream::commsTypes
            ) const override;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            ) override;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                scalarField& result,
                const bool add,
                const scalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const override;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<Type>& result,
                const bool add,
                const Field<Type>& psiInternal,
                const scalarField& coeffs,
                const Pstream::commsTypes commsType
            ) const override;


        // Block-coupled interface functionality

            //- Multiply the field by coefficients and add into the result
            void addToInternalCoupledField
            (
                Field<Type>&,
                const Field<Type>&,
                const CoeffField<Type>&
            ) const
            {
                NotImplemented;
            }

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const Field<Type>&,
                Field<Type>&,
                const BlockLduMatrix<Type>&,
                const CoeffField<Type>&,
                const Pstream::commsTypes commsType
            ) const override
            {
                NotImplemented;
            }


        // I-O

            //- Write
            virtual void write(Ostream& os) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fields/fvPatchFields/constraint/nonConformalDiscreteMixing/nonConformalDiscreteMixingFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
