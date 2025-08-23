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
    (c) 2011-2023 OpenFOAM Foundation
    (c) 2017-2025 Engys Ltd.

Class
    Foam::fvsPatchField

Description
    An abstract base class with a fat-interface to all derived classes
    covering all possible ways in which they might be used.

    The first level of derivation is to basic patchFields which cover
    zero-gradient, fixed-gradient, fixed-value and mixed conditions.

    The next level of derivation covers all the specialised typed with
    specific evaluation proceedures, particularly with respect to specific
    fields.

SourceFiles
    fvsPatchField.C
    fvsPatchFieldNew.C

\*---------------------------------------------------------------------------*/

#ifndef fvsPatchField_H
#define fvsPatchField_H

#include "fvMesh/fvPatches/fvPatch/fvPatch.H"
#include "fields/DimensionedFields/DimensionedField/DimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class objectRegistry;
class dictionary;
class fvPatchFieldMapper;
class surfaceMesh;
class gibFvPatchFieldMapper;


// Forward declaration of friend functions and operators

template<class Type>
class fvsPatchField;

template<class Type>
class calculatedFvsPatchField;

template<class Type>
Ostream& operator<<(Ostream&, const fvsPatchField<Type>&);


/*---------------------------------------------------------------------------*\
                        Class fvsPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fvsPatchField
:
    public Field<Type>
{
    // Private data

        //- Reference to patch
        const fvPatch& patch_;

        //- Reference to internal field
        const DimensionedField<Type, surfaceMesh>& internalField_;

        //- Optional patch type, used to allow specified boundary conditions
        //  to be applied to constraint patches by providing the constraint
        //  patch type as 'patchType'
        word patchType_;

        mutable autoPtr<Field<Type>> oldTimeFieldPtr_;

        Type defaultGIBValue_;

public:

    typedef fvPatch Patch;
    typedef calculatedFvsPatchField<Type> Calculated;


    //- Runtime type information
    TypeName("fvsPatchField");

    //- Debug switch to disallow the use of genericFvsPatchField
    static int disallowGenericFvsPatchField;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            fvsPatchField,
            patch,
            (
                const fvPatch& p,
                const DimensionedField<Type, surfaceMesh>& iF
            ),
            (p, iF)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            fvsPatchField,
            patchMapper,
            (
                const fvsPatchField<Type>& ptf,
                const fvPatch& p,
                const DimensionedField<Type, surfaceMesh>& iF,
                const fvPatchFieldMapper& m
            ),
            (dynamic_cast<const fvsPatchFieldType&>(ptf), p, iF, m)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            fvsPatchField,
            dictionary,
            (
                const fvPatch& p,
                const DimensionedField<Type, surfaceMesh>& iF,
                const dictionary& dict
            ),
            (p, iF, dict)
        );


    // Constructors

        //- Construct from patch and internal field
        fvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch and internal field and patch field
        fvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        fvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&,
            const bool valueRequired = true
        );

        //- Construct by mapping the given fvsPatchField onto a new patch
        fvsPatchField
        (
            const fvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired = true
        );

        //- Construct as copy
        fvsPatchField(const fvsPatchField<Type>&);

        //- Construct and return a clone
        virtual tmp<fvsPatchField<Type>> clone() const
        {
            return tmp<fvsPatchField<Type>>(new fvsPatchField<Type>(*this));
        }

        //- Construct as copy setting internal field reference
        fvsPatchField
        (
            const fvsPatchField<Type>&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvsPatchField<Type>> clone
        (
            const DimensionedField<Type, surfaceMesh>& iF
        ) const
        {
            return tmp<fvsPatchField<Type>>
            (
                new fvsPatchField<Type>(*this, iF)
            );
        }

        //- Construct and return a clone onto a new patch
        //  setting internal field reference
        virtual tmp<fvsPatchField<Type>> clone
        (
            const fvPatch& patch,
            const DimensionedField<Type, surfaceMesh>& iF
        ) const
        {
            return tmp<fvsPatchField<Type>>
            (
                new fvsPatchField<Type>(patch, iF, *this)
            );
        }


    // Selectors

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values)
        static tmp<fvsPatchField<Type>> New
        (
            const word&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values)
        //  Allows override of constraint type
        static tmp<fvsPatchField<Type>> New
        (
            const word&,
            const word& actualPatchType,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Return a pointer to a new patchField created on freestore from
        //  a given fvsPatchField mapped onto a new patch
        static tmp<fvsPatchField<Type>> New
        (
            const fvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Return a pointer to a new patchField created on freestore
        //  from dictionary
        static tmp<fvsPatchField<Type>> New
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Return a pointer to a new calculatedFvsPatchField created on
        //  freestore without setting patchField values
        static tmp<fvsPatchField<Type>> NewCalculatedType
        (
            const fvPatch&
        );

        //- Return a pointer to a new calculatedFvsPatchField created on
        //  freestore without setting patchField values
        template<class Type2>
        static tmp<fvsPatchField<Type>> NewCalculatedType
        (
            const fvsPatchField<Type2>&
        );


    //- Destructor
    virtual ~fvsPatchField() = default;


    // Member functions

        // Access

            //- Return local objectRegistry
            const objectRegistry& db() const;

            //- Return patch
            const fvPatch& patch() const
            {
                return patch_;
            }

            //- Return dimensioned internal field reference
            const DimensionedField<Type, surfaceMesh>& internalField() const
            {
                return internalField_;
            }

            //- Return internal field reference
            const Field<Type>& primitiveField() const
            {
                return internalField_;
            }

            //- Return old field reference
            const Field<Type>& oldTimeField() const
            {
                return oldTimeFieldPtr_();
            }

            const Type& defaultGIBValue() const
            {
                return defaultGIBValue_;
            }

            //- Optional patch type
            const word& patchType() const
            {
                return patchType_;
            }

            //- Optional patch type
            word& patchType()
            {
                return patchType_;
            }

            virtual word constraintType()
            {
                return word("");
            }

            //- Return the type of the calculated for of fvsPatchField
            static const word& calculatedType();

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return false;
            }

            //- Return true if this patch field is coupled
            virtual bool coupled() const
            {
                return false;
            }


        // Evaluation functions

            //- Initialise return of the patchField on the opposite patch of a
            //  coupled patch
            virtual void initPatchNeighbourField
            (
                const Pstream::commsTypes commsType
            ) const
            {}

            //- Return patchField on the opposite patch of a coupled patch
            virtual tmp<Field<Type>> patchNeighbourField
            (
                const Pstream::commsTypes commsType
            ) const
            {
                NotImplemented;
                return *this;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvsPatchField onto this fvsPatchField
            virtual void rmap
            (
                const fvsPatchField<Type>&,
                const labelList&
            );

            virtual void updateGIB();

            virtual void autoMapGIB(const gibFvPatchFieldMapper&);

            //- store the previous time of GIB
            virtual void storeGIB();

            //- Reset the fvsPatchField to the given fvsPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvsPatchField<Type>&);

            //- Return patch value for explicit gradient evaluation
            virtual tmp<Field<Type>> gradientBoundaryValue() const
            {
                return *this;
            }

        //- Write
        virtual void write(Ostream&) const;


        // Check

            //- Check fvsPatchField<Type> against given fvsPatchField<Type>
            void check(const fvsPatchField<Type>&) const;


    // Member operators

        virtual void operator=(const UList<Type>&);

        virtual void operator=(const fvsPatchField<Type>&);
        virtual void operator+=(const fvsPatchField<Type>&);
        virtual void operator-=(const fvsPatchField<Type>&);
        virtual void operator*=(const fvsPatchField<scalar>&);
        virtual void operator/=(const fvsPatchField<scalar>&);

        virtual void operator+=(const Field<Type>&);
        virtual void operator-=(const Field<Type>&);

        virtual void operator*=(const Field<scalar>&);
        virtual void operator/=(const Field<scalar>&);

        virtual void operator=(const Type&);
        virtual void operator+=(const Type&);
        virtual void operator-=(const Type&);
        virtual void operator*=(const scalar);
        virtual void operator/=(const scalar);


        // Force an assignment irrespective of form of patch

        virtual void forceAssign(const fvsPatchField<Type>&);
        virtual void forceAssign(const Field<Type>&);
        virtual void forceAssign(const Type&);


    // Ostream operator

        friend Ostream& operator<< <Type>(Ostream&, const fvsPatchField<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fields/fvsPatchFields/fvsPatchField/fvsPatchField.C"
    #include "fields/fvsPatchFields/basic/calculated/calculatedFvsPatchField.H"
#endif


#define makeFvsPatchTypeFieldTypeName(type)                                    \
    defineNamedTemplateTypeNameAndDebug(type, 0);


#define makeFvsPatchFieldsTypeName(type)                                       \
    makeFvsPatchTypeFieldTypeName(type##FvsPatchScalarField);                  \
    makeFvsPatchTypeFieldTypeName(type##FvsPatchVectorField);                  \
    makeFvsPatchTypeFieldTypeName(type##FvsPatchSphericalTensorField);         \
    makeFvsPatchTypeFieldTypeName(type##FvsPatchSymmTensorField);              \
    makeFvsPatchTypeFieldTypeName(type##FvsPatchTensorField);


// For templated patch fields
#define makeFvsPatchTypeField(PatchTypeField, typePatchTypeField)              \
    defineNamedTemplateTypeNameAndDebug(typePatchTypeField, 0);                \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        patch                                                                  \
    );                                                                         \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        patchMapper                                                            \
    );                                                                         \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        dictionary                                                             \
    );


// For non-templated patch fields
#define makeUntemplatedFvsPatchTypeField(PatchTypeField, typePatchTypeField)   \
    defineTypeNameAndDebug(typePatchTypeField, 0);                             \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        patch                                                                  \
    );                                                                         \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        patchMapper                                                            \
    );                                                                         \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        dictionary                                                             \
    );


#define makeFvsPatchFields(type)                                               \
    makeFvsPatchTypeField(fvsPatchScalarField, type##FvsPatchScalarField);     \
    makeFvsPatchTypeField(fvsPatchVectorField, type##FvsPatchVectorField);     \
    makeFvsPatchTypeField                                                      \
    (                                                                          \
        fvsPatchSphericalTensorField,                                          \
        type##FvsPatchSphericalTensorField                                     \
    );                                                                         \
    makeFvsPatchTypeField                                                      \
    (                                                                          \
        fvsPatchSymmTensorField,                                               \
        type##FvsPatchSymmTensorField                                          \
    );                                                                         \
    makeFvsPatchTypeField(fvsPatchTensorField, type##FvsPatchTensorField);


#define makeFvsPatchTypeFieldTypedef(fieldType, type)                          \
    typedef type##FvsPatchField<fieldType>                                     \
        CAT4(type, FvsPatch, CAPITALIZE(fieldType), Field);


#define makeFvsPatchTypeFieldTypedefs(type)                                    \
    typedef type##FvsPatchField<scalar> type##FvsPatchScalarField;             \
    typedef type##FvsPatchField<vector> type##FvsPatchVectorField;             \
    typedef type##FvsPatchField<sphericalTensor>                               \
        type##FvsPatchSphericalTensorField;                                    \
    typedef type##FvsPatchField<symmTensor> type##FvsPatchSymmTensorField;     \
    typedef type##FvsPatchField<tensor> type##FvsPatchTensorField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
