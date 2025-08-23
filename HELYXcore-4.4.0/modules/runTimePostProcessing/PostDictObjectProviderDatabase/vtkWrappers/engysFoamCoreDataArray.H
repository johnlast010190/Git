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
    (c) 2025 Engys Ltd.

Class
    engysFoamCoreDataArray

Description
    Base wrapper for arrays belonging to the openFOAM core

SourceFiles
    engysFoamCoreDataArray.txx

\*---------------------------------------------------------------------------*/

#pragma once

#include "vtkGenericDataArray.h"
#include "fields/Fields/fieldTypes.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
                   Class engysFoamCoreDataArray Declaration
\*---------------------------------------------------------------------------*/

template <class CoreValueTypeT, typename BaseValueTypeT>
class engysFoamCoreDataArray : public vtkGenericDataArray<engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>, BaseValueTypeT>
{
    typedef vtkGenericDataArray<engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>, BaseValueTypeT> GenericDataArrayType;

public:
    typedef engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT> SelfType;
    vtkAbstractTypeMacroWithNewInstanceType(SelfType, GenericDataArrayType, vtkDataArray, typeid(SelfType).name());

    typedef typename Superclass::ValueType ValueType;

    static engysFoamCoreDataArray* New();

    void SetCorePointer(const Foam::UList<CoreValueTypeT>* corePointer)
    {
        this->CorePointer = corePointer;
        this->CoreTupleSize = corePointer->size();
        this->CoreValueSize = corePointer->size() * this->CoreNumberOfComponents;
        this->Size = this->CoreValueSize;
        this->MaxId = this->Size - 1;
        this->Superclass::SetNumberOfComponents(this->CoreNumberOfComponents);
    }

    Foam::label GetCoreNumberOfTuples() const
    {
        return this->CoreTupleSize;
    }

    Foam::label GetCoreNumberOfValues() const
    {
        return this->CoreValueSize;
    }

    void SetNumberOfComponents(int nComponents) override
    {
        vtkErrorMacro("Cannot change the number of components of this class");
    }

    /**
   * Get the value at @a valueIdx
   */
    ValueType GetValue(vtkIdType valueIdx) const
    VTK_EXPECTS(0 <= valueIdx && valueIdx < this->GetNumberOfValues())
    {
        vtkIdType tupleIdx = valueIdx / this->CoreNumberOfComponents;
        int c = valueIdx % this->CoreNumberOfComponents;
        CoreValueTypeT coreValue = this->GetCoreValue(tupleIdx);
        return Foam::component(coreValue, c);
    }

    /**
   * Set the value at @a valueIdx to @a value
   */
    void SetValue(vtkIdType valueIdx, ValueType value)
    VTK_EXPECTS(this->GetCoreNumberOfValues() <= valueIdx && valueIdx < this->GetNumberOfValues())
    {
        assert(valueIdx >= this->GetCoreNumberOfValues() && valueIdx < this->GetNumberOfValues());
        vtkIdType tupleId = valueIdx / this->CoreNumberOfComponents;
        int c = valueIdx % this->CoreNumberOfComponents;
        Foam::setComponent((*this->ExtraBuffer)[tupleId - this->CoreTupleSize], c) = value;
    }

    ///@{
    /**
     * Copy the tuple at @a tupleIdx into @a tuple.
     */
    void GetTypedTuple(vtkIdType tupleIdx, ValueType* tuple) const
    VTK_EXPECTS(0 <= tupleIdx && tupleIdx < this->GetNumberOfTuples())
    {
        CoreValueTypeT coreValue = this->GetCoreValue(tupleIdx);
        for (Foam::direction c = 0; c < this->CoreNumberOfComponents; c++)
        {
            tuple[c] = Foam::component(coreValue, c);
        }
    }
    ///@}

    ///@{
    /**
     * Set this array's tuple at @a tupleIdx to the values in @a tuple.
     */
    void SetTypedTuple(vtkIdType tupleIdx, const ValueType* tuple)
    VTK_EXPECTS(this->GetCoreNumberOfTuples() <= tupleIdx && tupleIdx < this->GetNumberOfTuples())
    {
        assert(tupleIdx >= this->GetCoreNumberOfTuples() && tupleIdx < this->GetNumberOfTuples());
        CoreValueTypeT& coreValue = (*this->ExtraBuffer)[tupleIdx - this->CoreTupleSize];
        for (Foam::direction c = 0; c < this->CoreNumberOfComponents; c++)
        {
            Foam::setComponent(coreValue, c) = tuple[c];
        }
    }
    ///@}

    /**
     * Get component @a comp of the tuple at @a tupleIdx.
     */
    ValueType GetTypedComponent(vtkIdType tupleIdx, int comp) const
    VTK_EXPECTS(0 <= tupleIdx && tupleIdx < this->GetNumberOfTuples()) VTK_EXPECTS(0 <= comp && comp < this->GetNumberOfComponents())
    {
        CoreValueTypeT coreValue = this->GetCoreValue(tupleIdx);
        return Foam::component(coreValue, comp);
    }

    ///@{
    /**
     * Set component @a comp of the tuple at @a tupleIdx to @a value.
     */
    void SetTypedComponent(vtkIdType tupleIdx, int comp, ValueType value)
    VTK_EXPECTS(0 <= tupleIdx && tupleIdx < this->GetNumberOfTuples()) VTK_EXPECTS(0 <= comp && comp < this->GetNumberOfComponents())
    {
        assert(tupleIdx >= this->GetCoreNumberOfTuples() && tupleIdx < this->GetNumberOfTuples());
        CoreValueTypeT& coreValue = (*this->ExtraBuffer)[tupleIdx - this->CoreTupleSize];
        Foam::setComponent(coreValue, comp) = value;
    }
    ///@}

    ///@{
    /**
     * Set component @a comp of all extra writable tuples to @a value.
     */
    void FillTypedComponent(int compIdx, ValueType value) override;
    ///@}

    ///@{
    /**
     * Set all the values in the extra writable array to @a value.
     */
    void FillValue(ValueType value) override;
    void Fill(double value) override;
    ///@}

    // Overridden for optimized implementations:
    void SetTuple(vtkIdType tupleIdx, const float* tuple) override;
    void SetTuple(vtkIdType tupleIdx, const double* tuple) override;
    // MSVC doesn't like 'using' here (error C2487). Just forward instead:
    // using Superclass::SetTuple;
    void SetTuple(vtkIdType dstTupleIdx, vtkIdType srcTupleIdx, vtkAbstractArray* source) override
    {
        this->Superclass::SetTuple(dstTupleIdx, srcTupleIdx, source);
    }

    void GetTuple(vtkIdType tupleIdx, double* tuple) override;
    double* GetTuple(vtkIdType tupleIdx) override;

    /**
     * Tell the array explicitly that a single data element has
     * changed. Like DataChanged(), then is only necessary when you
     * modify the array contents without using the array's API.
     * @note This is a legacy method from vtkDataArrayTemplate, and is only
     * implemented for array-of-struct arrays. It currently just calls
     * DataChanged() and does nothing clever.
     * TODO this is only defined for AOS (vtkDataArrayTemplate leftover).
     * Deprecate to favor DataChanged?
     */
    void DataElementChanged(vtkIdType) { this->DataChanged(); }

    void SetNumberOfTuples(vtkIdType nTuples) override;

    [[nodiscard]] int GetArrayType() const override { return vtkAbstractArray::MappedDataArray; }
    VTK_NEWINSTANCE vtkArrayIterator* NewIterator() override { return nullptr; }
    [[nodiscard]] bool HasStandardMemoryLayout() const override { return false; }
    void DeepCopy(vtkDataArray* other) override;
    void ShallowCopy(vtkDataArray* other) override;

    // Reimplemented for efficiency:
    void InsertTuples(
        vtkIdType dstStart, vtkIdType n, vtkIdType srcStart, vtkAbstractArray* source) override;
    // MSVC doesn't like 'using' here (error C2487). Just forward instead:
    // using Superclass::InsertTuples;
    void InsertTuples(vtkIdList* dstIds, vtkIdList* srcIds, vtkAbstractArray* source) override;
    void InsertTuplesStartingAt(vtkIdType dstStart, vtkIdList* srcIds, vtkAbstractArray* source) override;

    void* GetVoidPointer(vtkIdType pos) override;

protected:
    engysFoamCoreDataArray();
    ~engysFoamCoreDataArray() override;

    VTK_NEWINSTANCE
    vtkObjectBase* NewInstanceInternal() const override { return vtkDataArray::CreateDataArray(Superclass::VTK_DATA_TYPE); }

    /**
     * Allocate space for numTuples. Old data is not preserved. If numTuples == 0,
     * all data is freed.
     */
    bool AllocateTuples(vtkIdType numTuples);

    /**
     * Allocate space for numTuples. Old data is preserved. If numTuples == 0,
     * all data is freed.
     */
    bool ReallocateTuples(vtkIdType numTuples);

    inline CoreValueTypeT GetCoreValue(vtkIdType tupleIdx) const;

    const Foam::UList<CoreValueTypeT>* CorePointer;
    std::shared_ptr<std::vector<CoreValueTypeT>> ExtraBuffer;
    std::shared_ptr<Foam::List<CoreValueTypeT>> VoidPointerBuffer;

    Foam::label CoreValueSize = 0;
    Foam::label CoreTupleSize = 0;
    static constexpr int CoreNumberOfComponents = Foam::pTraits<CoreValueTypeT>::nComponents;

    void InitializeVoidPointerVariable();

private:
    engysFoamCoreDataArray(const engysFoamCoreDataArray&) = delete;
    void operator=(const engysFoamCoreDataArray&) = delete;

    friend class vtkGenericDataArray<engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>, BaseValueTypeT>;
};


// ************************************************************************* //
