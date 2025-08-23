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
    (c) 2022-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#ifndef materialMacros_H
#define materialMacros_H

#include "cfdTools/general/include/fvCFD.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Define field functions
#define defineGeometricFieldFunction(Type,funcName)                           \
    virtual tmp<VolField<Type>> funcName##Geometric() const override          \
    {                                                                         \
        tmp<VolField<Type>> tField                                            \
        (                                                                     \
            VolField<Type>::New                                               \
            (                                                                 \
                #funcName,                                                    \
                obr_,                                                         \
                mesh(),                                                       \
                dimensioned<Type>                                             \
                (                                                             \
                    #funcName,                                                \
                    funcName##Model::dimensions(),                            \
                    Zero                                                      \
                )                                                             \
            )                                                                 \
        );                                                                    \
        VolField<Type>& field = tField.ref();                                 \
                                                                              \
        /* Calculate internal field */                                        \
        field.primitiveFieldRef() = funcName##Internal();                     \
                                                                              \
        /* Calculate bounary patches */                                       \
        forAll(field.boundaryField(), patchi)                                 \
        {                                                                     \
            field.boundaryFieldRef()[patchi].forceAssign                      \
            (                                                                 \
                funcName##Patch(patchi)                                       \
            );                                                                \
        }                                                                     \
        return tField;                                                        \
    }


#define materialsIsConstZero(funcName)                                        \
    virtual bool is##funcName##Const() const override                         \
    {                                                                         \
        return true;                                                          \
    }


#define materialsIsConstOne(funcName, field1Ptr)                              \
    virtual bool is##funcName##Const() const override                         \
    {                                                                         \
        return field1Ptr->isConst();                                          \
    }


#define materialsIsConstTwo(funcName, field1Ptr, field2Ptr)                   \
    virtual bool is##funcName##Const() const override                         \
    {                                                                         \
        return field1Ptr->isConst() && field2Ptr->isConst();                  \
    }


#define materialsIsConstThree(funcName, field1Ptr, field2Ptr, field3Ptr)      \
    virtual bool is##funcName##Const() const override                         \
    {                                                                         \
        return                                                                \
            field1Ptr->isConst()                                              \
         && field2Ptr->isConst()                                              \
         && field3Ptr->isConst();                                             \
    }


#define materialsIsConstFour(funcName, field1Ptr, field2Ptr, field3Ptr, field4Ptr)\
    virtual bool is##funcName##Const() const override                         \
    {                                                                         \
        return                                                                \
            field1Ptr->isConst()                                              \
         && field2Ptr->isConst()                                              \
         && field3Ptr->isConst()                                              \
         && field4Ptr->isConst();                                             \
    }


#define defineFieldFunctionsZero(Type,funcName)                               \
                                                                              \
    /* Calculate cell value */                                                \
    virtual Type funcName##Cell(const label celli) const override             \
    {                                                                         \
        return funcName();                                                    \
    }                                                                         \
                                                                              \
    /* Calculate value from set p and T */                                    \
    virtual Type funcName##Value(scalar p, scalar T) const override           \
    {                                                                         \
        return funcName();                                                    \
    }                                                                         \
                                                                              \
    /* Calculate internal field */                                            \
    virtual tmp<Type##Field> funcName##Internal() const override              \
    {                                                                         \
        return tmp<Type##Field>(new Type##Field(mesh().nCells(), funcName()));\
    }                                                                         \
                                                                              \
    /* Get patch fields references */                                         \
    virtual tmp<Type##Field> funcName##Patch                                  \
    (                                                                         \
        const label patchi                                                    \
    ) const override                                                          \
    {                                                                         \
        return                                                                \
            tmp<Type##Field>                                                  \
            (                                                                 \
                new Type##Field(mesh().boundary()[patchi].size(), funcName()) \
            );                                                                \
    }                                                                         \
                                                                              \
    defineGeometricFieldFunction(Type,funcName)


#define defineFieldFunctionsOne(Type,funcName,fieldPtr)                       \
                                                                              \
    /* Calculate cell value */                                                \
    virtual Type funcName##Cell(const label celli) const override             \
    {                                                                         \
        return funcName((*fieldPtr)[celli]);                                  \
    }                                                                         \
                                                                              \
    /* Calculate value from set p and T */                                    \
    virtual Type funcName##Value(scalar p, scalar T) const override           \
    {                                                                         \
        return funcName((*fieldPtr).value(p, T));                             \
    }                                                                         \
                                                                              \
    /* Calculate internal field */                                            \
    virtual tmp<Type##Field> funcName##Internal() const override              \
    {                                                                         \
        tmp<scalarField> tfield(fieldPtr->primitiveField());                  \
        scalarField& field = tfield.ref();                                    \
        tmp<Type##Field> tInternalField(new Type##Field(field.size()));       \
        Type##Field& fieldCells = tInternalField.ref();                       \
        forAll(fieldCells, celli)                                             \
        {                                                                     \
            fieldCells[celli] = funcName(field[celli]);                       \
        }                                                                     \
        return tInternalField;                                                \
    }                                                                         \
                                                                              \
    /* Get patch fields references */                                         \
    virtual tmp<Type##Field> funcName##Patch                                  \
    (                                                                         \
        const label patchi                                                    \
    ) const override                                                          \
    {                                                                         \
        tmp<scalarField> tfield(fieldPtr->boundaryField()[patchi]);           \
        scalarField& field = tfield.ref();                                    \
        tmp<Type##Field> tPatchField(new Type##Field(field.size()));          \
        Type##Field& patchField = tPatchField.ref();                          \
        forAll(patchField, facei)                                             \
        {                                                                     \
            patchField[facei] = funcName(field[facei]);                       \
        }                                                                     \
        return tPatchField;                                                   \
    }                                                                         \
                                                                              \
    defineGeometricFieldFunction(Type,funcName)


#define defineFieldFunctionsTwo(Type,funcName,field1Ptr,field2Ptr)            \
                                                                              \
    /* Calculate cell value */                                                \
    virtual Type funcName##Cell(const label celli) const override             \
    {                                                                         \
        return funcName((*field1Ptr)[celli], (*field2Ptr)[celli]);            \
    }                                                                         \
                                                                              \
    /* Calculate value from set p and T */                                    \
    virtual Type funcName##Value(scalar p, scalar T) const override           \
    {                                                                         \
        return funcName((*field1Ptr).value(p, T), (*field2Ptr).value(p, T));  \
    }                                                                         \
                                                                              \
    /* Calculate internal field */                                            \
    virtual tmp<Type##Field> funcName##Internal() const override              \
    {                                                                         \
        tmp<scalarField> tfield1(field1Ptr->primitiveField());                \
        tmp<scalarField> tfield2(field2Ptr->primitiveField());                \
        scalarField& field1 = tfield1.ref();                                  \
        scalarField& field2 = tfield2.ref();                                  \
        tmp<Type##Field> tInternalField(new Type##Field(mesh().nCells()));    \
        Type##Field& fieldCells = tInternalField.ref();                       \
        forAll(fieldCells, celli)                                             \
        {                                                                     \
            fieldCells[celli] = funcName(field1[celli], field2[celli]);       \
        }                                                                     \
        return tInternalField;                                                \
    }                                                                         \
                                                                              \
    /* Get patch fields references */                                         \
    virtual tmp<Type##Field> funcName##Patch                                  \
    (                                                                         \
        const label patchi                                                    \
    ) const override                                                          \
    {                                                                         \
        tmp<scalarField> tfield1(field1Ptr->boundaryField()[patchi]);         \
        tmp<scalarField> tfield2(field2Ptr->boundaryField()[patchi]);         \
        scalarField& field1 = tfield1.ref();                                  \
        scalarField& field2 = tfield2.ref();                                  \
        tmp<Type##Field> tPatchField(new Type##Field(field1.size()));         \
        Type##Field& patchField = tPatchField.ref();                          \
        forAll(patchField, facei)                                             \
        {                                                                     \
            patchField[facei] = funcName(field1[facei], field2[facei]);       \
        }                                                                     \
        return tPatchField;                                                   \
    }                                                                         \
                                                                              \
    defineGeometricFieldFunction(Type,funcName)


#define defineFieldFunctionsThree(Type,funcName,field1Ptr,field2Ptr,field3Ptr)\
                                                                              \
    /* Calculate cell value */                                                \
    virtual Type funcName##Cell(const label celli) const override             \
    {                                                                         \
        return                                                                \
            funcName                                                          \
            (                                                                 \
                (*field1Ptr)[celli],                                          \
                (*field2Ptr)[celli],                                          \
                (*field3Ptr)[celli]                                           \
            );                                                                \
    }                                                                         \
                                                                              \
    /* Calculate value from set p and T */                                    \
    virtual Type funcName##Value(scalar p, scalar T) const override           \
    {                                                                         \
        return                                                                \
            funcName                                                          \
            (                                                                 \
                (*field1Ptr).value(p, T),                                     \
                (*field2Ptr).value(p, T),                                     \
                (*field3Ptr).value(p, T)                                      \
            );                                                                \
    }                                                                         \
                                                                              \
    /* Calculate internal field */                                            \
    virtual tmp<Type##Field> funcName##Internal() const override              \
    {                                                                         \
        tmp<scalarField> tfield1(field1Ptr->primitiveField());                \
        tmp<scalarField> tfield2(field2Ptr->primitiveField());                \
        tmp<scalarField> tfield3(field3Ptr->primitiveField());                \
        scalarField& field1 = tfield1.ref();                                  \
        scalarField& field2 = tfield2.ref();                                  \
        scalarField& field3 = tfield3.ref();                                  \
        tmp<Type##Field> tInternalField(new Type##Field(mesh().nCells()));    \
        Type##Field& fieldCells = tInternalField.ref();                       \
        forAll(fieldCells, celli)                                             \
        {                                                                     \
            fieldCells[celli] =                                               \
                funcName                                                      \
                (                                                             \
                    field1[celli],                                            \
                    field2[celli],                                            \
                    field3[celli]                                             \
                );                                                            \
        }                                                                     \
        return tInternalField;                                                \
    }                                                                         \
                                                                              \
    /* Get patch fields references */                                         \
    virtual tmp<Type##Field> funcName##Patch                                  \
    (                                                                         \
        const label patchi                                                    \
    ) const override                                                          \
    {                                                                         \
        tmp<scalarField> tfield1(field1Ptr->boundaryField()[patchi]);         \
        tmp<scalarField> tfield2(field2Ptr->boundaryField()[patchi]);         \
        tmp<scalarField> tfield3(field3Ptr->boundaryField()[patchi]);         \
        scalarField& field1 = tfield1.ref();                                  \
        scalarField& field2 = tfield2.ref();                                  \
        scalarField& field3 = tfield3.ref();                                  \
        tmp<Type##Field> tPatchField(new Type##Field(field1.size()));         \
        Type##Field& patchField = tPatchField.ref();                          \
        forAll(patchField, facei)                                             \
        {                                                                     \
            patchField[facei] =                                               \
                funcName                                                      \
                (                                                             \
                    field1[facei],                                            \
                    field2[facei],                                            \
                    field3[facei]                                             \
                );                                                            \
        }                                                                     \
        return tPatchField;                                                   \
    }                                                                         \
                                                                              \
    defineGeometricFieldFunction(Type,funcName)


#define defineFieldFunctionsFour(Type,funcName,field1Ptr,field2Ptr,field3Ptr,field4Ptr)\
                                                                              \
    virtual Type funcName##Cell(const label celli) const override             \
    {                                                                         \
        return                                                                \
            funcName                                                          \
            (                                                                 \
                (*field1Ptr)[celli],                                          \
                (*field2Ptr)[celli],                                          \
                (*field3Ptr)[celli],                                          \
                (*field4Ptr)[celli]                                           \
            );                                                                \
    }                                                                         \
                                                                              \
    /* Calculate value from set p and T */                                    \
    virtual Type funcName##Value(scalar p, scalar T) const override           \
    {                                                                         \
        return                                                                \
            funcName                                                          \
            (                                                                 \
                (*field1Ptr).value(p, T),                                     \
                (*field2Ptr).value(p, T),                                     \
                (*field3Ptr).value(p, T)                                      \
                (*field4Ptr).value(p, T)                                      \
            );                                                                \
    }                                                                         \
                                                                              \
    /* Calculate internal field */                                            \
    virtual tmp<Type##Field> funcName##Internal() const override              \
    {                                                                         \
        tmp<scalarField> tfield1(field1Ptr->primitiveField());                \
        tmp<scalarField> tfield2(field2Ptr->primitiveField());                \
        tmp<scalarField> tfield3(field3Ptr->primitiveField());                \
        tmp<scalarField> tfield4(field4Ptr->primitiveField());                \
        scalarField& field1 = tfield1.ref();                                  \
        scalarField& field2 = tfield2.ref();                                  \
        scalarField& field3 = tfield3.ref();                                  \
        scalarField& field4 = tfield4.ref();                                  \
        tmp<Type##Field> tInternalField(new Type##Field(mesh().nCells()));    \
        Type##Field& fieldCells = tInternalField.ref();                       \
        forAll(fieldCells, celli)                                             \
        {                                                                     \
            fieldCells[celli] =                                               \
                funcName                                                      \
                (                                                             \
                    field1[celli],                                            \
                    field2[celli],                                            \
                    field3[celli],                                            \
                    field4[celli]                                             \
                );                                                            \
        }                                                                     \
        return tInternalField;                                                \
    }                                                                         \
                                                                              \
    /* Get patch fields references */                                         \
    virtual tmp<Type##Field> funcName##Patch                                  \
    (                                                                         \
        const label patchi                                                    \
    ) const override                                                          \
    {                                                                         \
        tmp<scalarField> tfield1(field1Ptr->boundaryField()[patchi]);         \
        tmp<scalarField> tfield2(field2Ptr->boundaryField()[patchi]);         \
        tmp<scalarField> tfield3(field3Ptr->boundaryField()[patchi]);         \
        tmp<scalarField> tfield4(field4Ptr->boundaryField()[patchi]);         \
        scalarField& field1 = tfield1.ref();                                  \
        scalarField& field2 = tfield2.ref();                                  \
        scalarField& field3 = tfield3.ref();                                  \
        scalarField& field4 = tfield4.ref();                                  \
        tmp<Type##Field> tPatchField(new Type##Field(field1.size()));         \
        Type##Field& patchField = tPatchField.ref();                          \
        forAll(patchField, facei)                                             \
        {                                                                     \
            patchField[facei] =                                               \
                funcName                                                      \
                (                                                             \
                    field1[facei],                                            \
                    field2[facei],                                            \
                    field3[facei],                                            \
                    field4[facei]                                             \
                );                                                            \
        }                                                                     \
        return tPatchField;                                                   \
    }                                                                         \
                                                                              \
    defineGeometricFieldFunction(Type,funcName)


#define declareTypeFunctionZero(returnType, funcName)                         \
    inline returnType funcName() const;


#define declareTypeFunctionOne(returnType, funcName, Type1)                   \
    inline returnType funcName(const Type1 val) const;


#define declareTypeFunctionTwo(returnType, funcName, Type1, Type2)            \
    inline returnType funcName(const Type1 val1, const Type2 val2) const;


#define declareTypeFunctionThree(returnType,funcName,Type1,Type2,Type3)       \
    inline returnType funcName                                                \
    (                                                                         \
        const Type1 val1,                                                     \
        const Type2 val2,                                                     \
        const Type3 val3                                                      \
    ) const;


#define declareTypeFunctionFour(returnType,funcName,Type1,Type2,Type3,Type4)  \
    inline returnType funcName                                                \
    (                                                                         \
        const Type1 val1,                                                     \
        const Type2 val2,                                                     \
        const Type3 val3,                                                     \
        const Type4 val4                                                      \
    ) const;


#define scalarFieldFunctionsNone(funcName)                                    \
    materialsIsConstZero(funcName)                                            \
    inline scalar funcName() const { return 0.0; }                            \
    defineFieldFunctionsZero(scalar,funcName)


#define vectorFieldFunctionsNone(funcName)                                    \
    materialsIsConstZero(funcName)                                            \
    inline vector funcName() const { return vector::zero; }                   \
    defineFieldFunctionsZero(vector,funcName)


#define tensorFieldFunctionsNone(funcName)                                    \
    materialsIsConstZero(funcName)                                            \
    inline tensor funcName() const { return tensor::zero; }                   \
    defineFieldFunctionsZero(tensor,funcName)


#define scalarFieldFunctionsZero(funcName)                                    \
    materialsIsConstZero(funcName)                                            \
    declareTypeFunctionZero(scalar, funcName)                                 \
    defineFieldFunctionsZero(scalar,funcName)


#define vectorFieldFunctionsZero(funcName)                                    \
    materialsIsConstZero(funcName)                                            \
    declareTypeFunctionZero(vector, funcName)                                 \
    defineFieldFunctionsZero(vector,funcName)


#define tensorFieldFunctionsZero(funcName)                                    \
    materialsIsConstZero(funcName)                                            \
    declareTypeFunctionZero(tensor, funcName)                                 \
    defineFieldFunctionsZero(tensor,funcName)


#define scalarFieldFunctionsOne(funcName,fieldPtr)                            \
    materialsIsConstOne(funcName,fieldPtr)                                    \
    declareTypeFunctionOne(scalar,funcName,scalar)                            \
    defineFieldFunctionsOne(scalar,funcName,fieldPtr)


#define vectorFieldFunctionsOne(funcName,fieldPtr)                            \
    materialsIsConstOne(funcName,fieldPtr)                                    \
    declareTypeFunctionOne(vector,funcName,scalar)                            \
    defineFieldFunctionsOne(vector,funcName,fieldPtr)


#define tensorFieldFunctionsOne(funcName,fieldPtr)                            \
    materialsIsConstOne(funcName,fieldPtr)                                    \
    declareTypeFunctionOne(tensor, funcName, scalar)                          \
    defineFieldFunctionsOne(tensor,funcName,fieldPtr)


#define scalarFieldFunctionsTwo(funcName,field1Ptr,field2Ptr)                 \
    materialsIsConstTwo(funcName, field1Ptr, field2Ptr)                       \
    declareTypeFunctionTwo(scalar, funcName, scalar, scalar)                  \
    defineFieldFunctionsTwo(scalar,funcName,field1Ptr,field2Ptr)


#define vectorFieldFunctionsTwo(funcName,field1Ptr,field2Ptr)                 \
    materialsIsConstTwo(funcName, field1Ptr, field2Ptr)                       \
    declareTypeFunctionTwo(vector, funcName, scalar, scalar)                  \
    defineFieldFunctionsTwo(vector,funcName,field1Ptr,field2Ptr)


#define tensorFieldFunctionsTwo(funcName,field1Ptr,field2Ptr)                 \
    materialsIsConstTwo(funcName, field1Ptr, field2Ptr)                       \
    declareTypeFunctionTwo(tensor, funcName, scalar, scalar)                  \
    defineFieldFunctionsTwo(tensor,funcName,field1Ptr,field2Ptr)


#define scalarFieldFunctionsThree(funcName,field1Ptr,field2Ptr,field3Ptr)\
    materialsIsConstThree(funcName, field1Ptr, field2Ptr, field3Ptr)          \
    declareTypeFunctionThree(scalar, funcName, scalar, scalar, scalar)        \
    defineFieldFunctionsThree(scalar,funcName,field1Ptr,field2Ptr,field3Ptr)


#define vectorFieldFunctionsThree(funcName,field1Ptr,field2Ptr,field3Ptr)\
    materialsIsConstThree(funcName, field1Ptr, field2Ptr, field3Ptr)          \
    declareTypeFunctionThree(vector, funcName, scalar, scalar, scalar)        \
    defineFieldFunctionsThree(vector,funcName,field1Ptr,field2Ptr,field3Ptr)


#define tensorFieldFunctionsThree(funcName,field1Ptr,field2Ptr,field3Ptr)     \
    materialsIsConstThree(funcName, field1Ptr, field2Ptr, field3Ptr)          \
    declareTypeFunctionThree(tensor, funcName, scalar, scalar, scalar)        \
    defineFieldFunctionsThree(tensor,funcName,field1Ptr,field2Ptr,field3Ptr)


#define scalarFieldFunctionsFour(funcName,field1Ptr,field2Ptr,field3Ptr,field4Ptr)\
    materialsIsConstFour(funcName,field1Ptr,field2Ptr,field3Ptr,field4Ptr)    \
    declareTypeFunctionFour(scalar,funcName,scalar,scalar,scalar,scalar)      \
    defineFieldFunctionsFour(scalar,funcName,field1Ptr,field2Ptr,field3Ptr,field4Ptr)


#define vectorFieldFunctionsFour(funcName,field1Ptr,field2Ptr,field3Ptr,field4Ptr)\
    materialsIsConstFour(funcName,field1Ptr,field2Ptr,field3Ptr,field4Ptr)    \
    declareTypeFunctionFour(vector,funcName,scalar,scalar,scalar,scalar)      \
    defineFieldFunctionsFour(vector,funcName,field1Ptr,field2Ptr,field3Ptr,field4Ptr)


#define tensorFieldFunctionsFour(funcName,field1Ptr,field2Ptr,field3Ptr,field4Ptr)\
    materialsIsConstFour(funcName,field1Ptr,field2Ptr,field3Ptr,field4Ptr)    \
    declareTypeFunctionFour(tensor,funcName,scalar,scalar,scalar,scalar)      \
    defineFieldFunctionsFour(tensor,funcName,field1Ptr,field2Ptr,field3Ptr,field4Ptr)


// Funcion definitioins no const functions
#define scalarFieldFunctionsNoConstNone(funcName)                             \
    inline scalar funcName() const { return 0.0; }                            \
    defineFieldFunctionsZero(scalar,funcName)


#define scalarFieldFunctionsNoConstZero(funcName)                             \
    declareTypeFunctionZero(scalar, funcName)                                 \
    defineFieldFunctionsZero(scalar, funcName)


#define scalarFieldFunctionsNoConstOne(funcName,field1Ptr)                    \
    declareTypeFunctionOne(scalar, funcName, scalar)                          \
    defineFieldFunctionsOne(scalar, funcName, field1Ptr)


#define scalarFieldFunctionsNoConstTwo(funcName,field1Ptr,field2Ptr)          \
    declareTypeFunctionTwo(scalar, funcName, scalar, scalar)                  \
    defineFieldFunctionsTwo(scalar, funcName, field1Ptr, field2Ptr)


#define scalarFieldFunctionsNoConstThree(funcName,field1Ptr,field2Ptr,field3Ptr)\
    declareTypeFunctionThree(scalar,funcName,scalar,scalar,scalar)            \
    defineFieldFunctionsThree(scalar,funcName,field1Ptr,field2Ptr,field3Ptr)


#define scalarFieldFunctionsNoConstFour(funcName,field1Ptr,field2Ptr,field3Ptr,field4Ptr)\
    declareTypeFunctionFour(scalar,funcName,scalar,scalar,scalar,scalar)      \
    defineFieldFunctionsFour(scalar,funcName,field1Ptr,field2Ptr,field3Ptr,field4Ptr)


//- Declare functions
#define declareVirtualFieldFunctions(Type,funcName,virtual)                   \
    virtual Type funcName##Cell(const label celli) const override;            \
                                                                              \
    virtual Type funcName##Value(scalar p, scalar T) const override;          \
                                                                              \
    virtual tmp<Type##Field> funcName##Internal() const override;             \
                                                                              \
    virtual tmp<Type##Field> funcName##Patch                                  \
    (                                                                         \
        const label patchi                                                    \
    ) const override;                                                         \
                                                                              \
    virtual tmp<VolField<Type>> funcName##Geometric() const override;


#define declareFieldFunctions(Type,funcName)                                  \
    declareVirtualFieldFunctions(Type,funcName,)


#define declareScalarFieldFunctions(funcName)                                 \
    declareVirtualFieldFunctions(scalar,funcName,)


#define declareVectorFieldFunctions(funcName)                                 \
    declareVirtualFieldFunctions(vector,funcName,)


#define declareTensorFieldFunctions(funcName)                                 \
    declareVirtualFieldFunctions(tensor,funcName,)


#define defineMaterialModel(funcName,type,PatchField,GeoMesh)                 \
    namespace Foam                                                            \
    {                                                                         \
    class funcName##Model                                                     \
    :                                                                         \
        public baseModels<type,PatchField,GeoMesh>                            \
    {                                                                         \
    public:                                                                   \
        TypeName(#funcName);                                                  \
        virtual ~funcName##Model(){delete boundary_;};                        \
        static const dimensionSet modelDims;                                  \
        virtual bool mixture() const override { return false; }               \
        virtual bool incompressible() const override {return false;}          \
        virtual bool isochoric() const override { return false; }             \
        virtual bool isotropic() const override { return true; }              \
        static const bool isMixture;                                          \
        class Boundary                                                        \
        :                                                                     \
            public baseModels<type,PatchField,GeoMesh>::Boundary              \
        {                                                                     \
        public:                                                               \
            Boundary(const funcName##Model& base):funcName##Model_(base)      \
            {};                                                               \
            const funcName##Model& funcName##Model_;                          \
            virtual ~Boundary(){};                                            \
            virtual tmp<Field<type>> operator[](const label patchi) const     \
            {                                                                 \
                return funcName##Model_.funcName##Patch(patchi);              \
            }                                                                 \
        };                                                                    \
        typedef type ModelType;                                               \
        static const word tempName;                                           \
        Boundary* boundary_;                                                  \
        funcName##Model(){boundary_ = new Boundary(*this);};                  \
        virtual bool is##funcName##Const() const = 0;                         \
        virtual tmp<GeometricField<type,PatchField,GeoMesh>>                  \
        funcName##Geometric() const = 0;                                      \
        virtual tmp<type##Field> funcName##Internal() const = 0;              \
        virtual tmp<type##Field> funcName##Patch                              \
        (                                                                     \
            const label patchi                                                \
        ) const = 0;                                                          \
        virtual type funcName##Cell(const label celli) const = 0;             \
        virtual type funcName##Value(scalar p, scalar T) const = 0;           \
        virtual type value(scalar p, scalar T) const override                 \
        {                                                                     \
            return funcName##Value(p, T);                                     \
        }                                                                     \
        virtual bool isConst() const override                                 \
        {                                                                     \
            return is##funcName##Const();                                     \
        }                                                                     \
        virtual tmp<type##Field> primitiveField() const override              \
        {                                                                     \
            return funcName##Internal();                                      \
        }                                                                     \
        virtual baseModels<type,PatchField,GeoMesh>::Boundary&                \
        boundaryField() const override                                        \
        {                                                                     \
            return *boundary_;                                                \
        }                                                                     \
        virtual const dimensionSet& dimensions() const override               \
        {                                                                     \
            return modelDims;                                                 \
        }                                                                     \
        virtual type operator[](const label celli) const override             \
        {                                                                     \
            return funcName##Cell(celli);                                     \
        }                                                                     \
        virtual tmp<GeometricField<type,PatchField,GeoMesh>>                  \
        operator()() const override                                           \
        {                                                                     \
            return funcName##Geometric();                                     \
        }                                                                     \
        virtual void updateTable(const word&) override = 0;                   \
        virtual const depList& dep() const override = 0;                      \
        virtual const word funcType() const override                          \
        {                                                                     \
            return #funcName;                                                 \
        }                                                                     \
        virtual const word modelType() const override                         \
        {                                                                     \
            return #type;                                                     \
        }                                                                     \
        virtual bool read() override = 0;                                     \
    };                                                                        \
    }


#define declareMaterialModelDims(funcName,dim,type,PatchField,GeoMesh,mix)    \
    defineTypeNameAndDebug(funcName##Model, 0);                               \
    const dimensionSet funcName##Model::modelDims = dim;                      \
    const bool funcName##Model::isMixture = mix;                              \
    baseModels<type,PatchField,GeoMesh>::addToTable<funcName##Model>          \
        test##funcName##_;                                                    \
    const word funcName##Model::tempName = #type;


#define defineScalarMaterialModel(funcName)                                   \
    defineMaterialModel(funcName,scalar,fvPatchField,volMesh)


#define defineVectorMaterialModel(funcName)                                   \
    defineMaterialModel(funcName,vector,fvPatchField,volMesh)


#define defineTensorMaterialModel(funcName)                                   \
    defineMaterialModel(funcName,tensor,fvPatchField,volMesh)


#define declareMaterialModelScalarMixture(funcName,dims)                      \
    declareMaterialModelDims(funcName,dims,scalar,fvPatchField,volMesh,true)


#define declareMaterialModelVectorMixture(funcName,dims)                      \
    declareMaterialModelDims(funcName,dims,vector,fvPatchField,volMesh,true)


#define declareMaterialModelTensorMixture(funcName,dims)                      \
    declareMaterialModelDims(funcName,dims,tensor,fvPatchField,volMesh,true)


#define declareMaterialModelScalarDims(funcName,dims)                         \
    declareMaterialModelDims(funcName,dims,scalar,fvPatchField,volMesh,false)


#define declareMaterialModelVectorDims(funcName,dims)                         \
    declareMaterialModelDims(funcName,dims,vector,fvPatchField,volMesh,false)


#define declareMaterialModelTensorDims(funcName,dims)                         \
    declareMaterialModelDims(funcName,dims,tensor,fvPatchField,volMesh,false)


#define addMatThermo(BaseThermo, Mixture, CThermo, whichTable)                \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        whichTable,                                                           \
        CThermo##BaseThermo##Mixture,                                         \
        objectRegistry                                                        \
    )


#define makeMatThermo(BaseThermo, Mixture, CThermo, matTypeName)              \
    typedef CThermo<BaseThermo::composite, Mixture>                           \
        CThermo##BaseThermo##Mixture;                                         \
                                                                              \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        CThermo##BaseThermo##Mixture,                                         \
        #matTypeName,                                                         \
        0                                                                     \
    );                                                                        \
    addMatThermo(BaseThermo, Mixture, CThermo, basicThermo);                  \
    addMatThermo(BaseThermo, Mixture, CThermo, fluidThermo);                  \
    addMatThermo(BaseThermo, Mixture, CThermo, BaseThermo)


// Lazy approach without splitting multiphase thermo
#define makeMatThermoNoComposite(BaseThermo, Mixture, CThermo, matTypeName)    \
    typedef CThermo<BaseThermo, Mixture> CThermo##BaseThermo##Mixture;        \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        CThermo##BaseThermo##Mixture,                                         \
        #matTypeName,                                                         \
        0                                                                     \
    );                                                                        \
    addMatThermo(BaseThermo, Mixture, CThermo, basicThermo);                  \
    addMatThermo(BaseThermo, Mixture, CThermo, fluidThermo);                  \
    addMatThermo(BaseThermo, Mixture, CThermo, BaseThermo)


#define makeMatSolidThermo(BaseThermo, CThermo, Mixture)                      \
    typedef BasicThermo<BaseThermo::composite, Mixture> BasicThermo##Mixture; \
    typedef CThermo<BaseThermo::composite, Mixture>                           \
        CThermo##BaseThermo##Mixture;                                         \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        CThermo##BaseThermo##Mixture,                                         \
        "solid",                                                              \
        0                                                                     \
    );                                                                        \
    addMatThermo(BaseThermo, Mixture, CThermo, basicThermo);                  \
    addMatThermo(BaseThermo, Mixture, CThermo, BaseThermo);                   \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        BaseThermo,                                                           \
        CThermo##BaseThermo##Mixture,                                         \
        dictionary                                                            \
    );


#define castMaterial(modelName, funcName)                                     \
    if (modelName == funcName##Model::typeName)                               \
    {                                                                         \
        return dynamic_cast<funcName##Model*>(this);                          \
    }


#define materialDepList()                                                     \
    virtual const depList& dep() const override                               \
    {                                                                         \
        return materialModel::dep();                                          \
    }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //