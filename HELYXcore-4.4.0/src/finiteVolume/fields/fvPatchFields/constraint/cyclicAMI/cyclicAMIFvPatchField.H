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
    (c) 2010-2023 Engys Ltd.

Class
    Foam::cyclicAMIFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition enforces a cyclic condition between a pair of
    boundaries, whereby communication between the patches is performed using
    an arbitrary mesh interface (AMI) interpolation.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            cyclicAMI;
    }
    \endverbatim

Note
    The outer boundary of the patch pairs must be similar, i.e. if the owner
    patch is transformed to the neighbour patch, the outer perimiter of each
    patch should be identical (or very similar).

See also
    Foam::AMIInterpolation

SourceFiles
    cyclicAMIFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicAMIFvPatchField_H
#define cyclicAMIFvPatchField_H

#include "fields/fvPatchFields/basic/coupled/coupledFvPatchField.H"
#include "AMIInterpolation/patches/cyclicAMI/cyclicAMILduInterfaceField/cyclicAMILduInterfaceField.H"
#include "fvMesh/fvPatches/constraint/cyclicAMI/cyclicAMIFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class cyclicAMIFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cyclicAMIFvPatchField
:
    virtual public cyclicAMILduInterfaceField,
    public coupledFvPatchField<Type>
{
    // Private data

        //- Local reference cast into the cyclic patch
        const cyclicAMIFvPatch& cyclicAMIPatch_;

    // Private Member Functions

        //- Do implicit coeffs or not
        bool implicitCoeffs () const;

        //- Calculate and return patch neighbour field
        tmp<Field<Type>> interpPatchNeighbourField() const;

        //- Return neighbour side field given internal fields
        template<class Type2>
        tmp<Field<Type2>> neighbourSideField
        (
            const Field<Type2>&
        ) const;


public:

    //- Runtime type information
    TypeName(cyclicAMIFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        cyclicAMIFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        cyclicAMIFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given cyclicAMIFvPatchField onto a new patch
        cyclicAMIFvPatchField
        (
            const cyclicAMIFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        cyclicAMIFvPatchField(const cyclicAMIFvPatchField<Type>&);

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new cyclicAMIFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        cyclicAMIFvPatchField
        (
            const cyclicAMIFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new cyclicAMIFvPatchField<Type>(*this, iF)
            );
        }

        virtual ~cyclicAMIFvPatchField ()
        {};

    // Member functions

        // Access

            //- Return local reference cast into the cyclic AMI patch
            virtual const cyclicAMIFvPatch& cyclicAMIPatch() const
            {
                return cyclicAMIPatch_;
            }

            virtual word constraintType()
            {
                return cyclicAMIFvPatch::typeName_();
            }


        // Evaluation functions

            //- Return true if coupled. Note that the underlying patch
            //  is not coupled() - the points don't align.
            virtual bool coupled() const;

            //- Return neighbour coupled internal cell data
            virtual tmp<Field<Type>> patchNeighbourField
            (
                const Pstream::commsTypes
            ) const;

            //- Return internal field next to patch as patch field
            virtual tmp<Field<Type>> patchInternalField() const;

            //- Return reference to neighbour patchField
            const cyclicAMIFvPatchField<Type>& nbrPatchField() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                scalarField& result,
                const bool add,
                const scalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<Type>&,
                const bool add,
                const Field<Type>&,
                const scalarField&,
                const Pstream::commsTypes commsType
            ) const;


        // Cyclic AMI coupled interface functions

            //- Return transformation between the coupled patches
            virtual const transformer& transform() const
            {
                return cyclicAMIPatch_.transform();
            }

            //- Return rank of component for transform
            virtual int rank() const
            {
                return pTraits<Type>::rank;
            }

        // Block coupled interface functionality


            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const Field<Type>&,
                Field<Type>&,
                const BlockLduMatrix<Type>&,
                const CoeffField<Type>&,
                const Pstream::commsTypes commsType
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
            ) const
            {
                NotImplemented;
            }

            virtual void addToInternalCoupledField
            (
                Field<Type>&,
                const Field<Type>&,
                const CoeffField<Type>&
            ) const
            {
                NotImplemented;
            }

            virtual void addToInternalCoupledCoveredField
            (
                Field<Type>&,
                const Field<Type>&,
                const CoeffField<Type>&
            ) const
            {
                NotImplemented;
            }

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientInternalCoeffs
            (
                const scalarField& deltaCoeffs
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientBoundaryCoeffs
            (
                const scalarField& deltaCoeffs
            ) const;

        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
