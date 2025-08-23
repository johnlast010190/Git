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

\*---------------------------------------------------------------------------*/

#pragma once
#include "engysFoamCoreDataArray.H"

#include "vtkObjectFactory.h"

template<class CoreValueTypeT, typename BaseValueTypeT>
engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT> *engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>::New()
{
    auto result = new engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>();
    result->InitializeObjectBase();
    return result;
}

template<class CoreValueTypeT, typename BaseValueTypeT>
engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>::engysFoamCoreDataArray()
{
    this->CorePointer = nullptr;
    this->ExtraBuffer = std::shared_ptr<std::vector<CoreValueTypeT>>(new std::vector<CoreValueTypeT>);
    this->CoreValueSize = 0;
    this->CoreTupleSize = 0;
    this->NumberOfComponents = this->CoreNumberOfComponents;
}

template<class CoreValueTypeT, typename BaseValueTypeT>
engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>::~engysFoamCoreDataArray()
{
    this->CorePointer = nullptr;
}

template<class CoreValueTypeT, typename BaseValueTypeT>
CoreValueTypeT engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>::GetCoreValue(vtkIdType tupleIdx) const
{
    assert(tupleIdx >= 0 && tupleIdx < this->GetNumberOfTuples());
    return (tupleIdx < this->CoreTupleSize) ?
           (*this->CorePointer)[tupleIdx] :
           (*this->ExtraBuffer)[tupleIdx - this->CoreTupleSize];
}

template<class CoreValueTypeT, typename BaseValueTypeT>
void engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>::FillTypedComponent(int compIdx, ValueType value)
{
    for (CoreValueTypeT& coreValue : (*this->ExtraBuffer))
    {
        Foam::setComponent(coreValue, compIdx) = value;
    }
    this->Modified();
}

template<class CoreValueTypeT, typename BaseValueTypeT>
void engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>::FillValue(ValueType value)
{
    for (CoreValueTypeT& coreValue : (*this->ExtraBuffer))
    {
        for (Foam::direction c = 0; c < this->CoreNumberOfComponents; c++)
        {
            Foam::setComponent(coreValue, c) = value;
        }
    }
    this->Modified();
}

template<class CoreValueTypeT, typename BaseValueTypeT>
void engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>::Fill(double value)
{
    for (CoreValueTypeT& coreValue : (*this->ExtraBuffer))
    {
        for (Foam::direction c = 0; c < this->CoreNumberOfComponents; c++)
        {
            Foam::setComponent(coreValue, c) = value;
        }
    }
    this->Modified();
}

template<class CoreValueTypeT, typename BaseValueTypeT>
void engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>::SetTuple(vtkIdType tupleIdx, const float *tuple)
{
    CoreValueTypeT& coreValue = (*this->ExtraBuffer)[tupleIdx - this->CoreTupleSize];
    for (Foam::direction c = 0; c < this->CoreNumberOfComponents; c++)
    {
        Foam::setComponent(coreValue, c) = tuple[c];
    }
    this->Modified();
}

template<class CoreValueTypeT, typename BaseValueTypeT>
void engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>::SetTuple(vtkIdType tupleIdx, const double *tuple)
{
    CoreValueTypeT& coreValue = (*this->ExtraBuffer)[tupleIdx - this->CoreTupleSize];
    for (Foam::direction c = 0; c < this->CoreNumberOfComponents; c++)
    {
        Foam::setComponent(coreValue, c) = tuple[c];
    }
    this->Modified();
}

//template<class CoreValueTypeT, typename BaseValueTypeT>
//void engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>::InsertTuple(vtkIdType tupleIdx, const float *tuple)
//{
//    vtkIdType extraTupleIdx = tupleIdx - CoreTupleSize;
//    if (extraTupleIdx >= (*this->ExtraBuffer).size()) (*this->ExtraBuffer).resize(extraTupleIdx + 1);
//    this->SetTuple(tupleIdx, tuple);
//}
//
//template<class CoreValueTypeT, typename BaseValueTypeT>
//void engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>::InsertTuple(vtkIdType tupleIdx, const double *tuple)
//{
//    vtkIdType extraTupleIdx = tupleIdx - CoreTupleSize;
//    if (extraTupleIdx >= this->ExtraBuffer->size()) this->ExtraBuffer->resize(extraTupleIdx + 1);
//    this->SetTuple(tupleIdx, tuple);
//}
//
//template<class CoreValueTypeT, typename BaseValueTypeT>
//void engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>::InsertComponent(vtkIdType tupleIdx, int compIdx, double value)
//{
//    vtkIdType extraTupleIdx = tupleIdx - CoreTupleSize;
//    if (extraTupleIdx >= this->ExtraBuffer->size()) this->ExtraBuffer->resize(extraTupleIdx + 1);
//    this->SetComponent(tupleIdx, compIdx, value);
//}

template<class CoreValueTypeT, typename BaseValueTypeT>
bool engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>::AllocateTuples(vtkIdType numTuples)
{
    if (numTuples >= this->CoreTupleSize)
    {
        size_t newExtraSize = numTuples - this->CoreTupleSize;
        this->ExtraBuffer->resize(newExtraSize);
        this->Size = numTuples * this->CoreNumberOfComponents;
        // Update MaxId if we truncated:
        if ((this->Size - 1) < this->MaxId)
        {
            this->MaxId = (this->Size - 1);
        }
        return true;
    }
    else
    {
        return false;
    }
}

template<class CoreValueTypeT, typename BaseValueTypeT>
bool engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>::ReallocateTuples(vtkIdType numTuples)
{
    if (numTuples >= this->CoreTupleSize)
    {
        size_t newExtraSize = numTuples - this->CoreTupleSize;
        this->ExtraBuffer->resize(newExtraSize);
        this->Size = numTuples * this->CoreNumberOfComponents;
        // Update MaxId if we truncated:
        if ((this->Size - 1) < this->MaxId)
        {
            this->MaxId = (this->Size - 1);
        }
        return true;
    }
    else
    {
        return false;
    }
}

template<class CoreValueTypeT, typename BaseValueTypeT>
double* engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>::GetTuple(vtkIdType tupleIdx)
{
    CoreValueTypeT coreValue = this->GetCoreValue(tupleIdx);
    for (Foam::direction c = 0; c < this->CoreNumberOfComponents; c++)
    {
        this->LegacyTuple[c] = Foam::component(coreValue, c);
    }
    return &this->LegacyTuple[0];
}

template<class CoreValueTypeT, typename BaseValueTypeT>
void engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>::GetTuple(vtkIdType tupleIdx, double* tuple)
{
    CoreValueTypeT coreValue = this->GetCoreValue(tupleIdx);
    for (Foam::direction c = 0; c < this->CoreNumberOfComponents; c++)
    {
        tuple[c] = Foam::component(coreValue, c);
    }
}

template<class CoreValueTypeT, typename BaseValueTypeT>
void engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>::DeepCopy(vtkDataArray* other)
{
    if (other == nullptr)
    {
        return;
    }

    if (this != other)
    {
        auto *o = dynamic_cast<engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT> *>(other);
        if (o)
        {
            this->CorePointer = o->CorePointer;
            this->CoreTupleSize = o->CoreTupleSize;
            this->CoreValueSize = o->CoreValueSize;
            if (o->VoidPointerBuffer)
            {
                this->VoidPointerBuffer = std::shared_ptr<Foam::List<CoreValueTypeT>>(new Foam::List<CoreValueTypeT>(*o->VoidPointerBuffer.get()));
                this->CorePointer = this->VoidPointerBuffer.get();
            }
            this->ExtraBuffer->assign(o->ExtraBuffer->begin(), o->ExtraBuffer->end());
            this->Size = o->Size;
            this->MaxId = o->MaxId;
            this->SetLookupTable(nullptr);
            this->Modified();
        }
    }
}

template<class CoreValueTypeT, typename BaseValueTypeT>
void engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>::ShallowCopy(vtkDataArray* other)
{
    if (other == nullptr)
    {
        return;
    }

    if (this != other)
    {
        auto *o = dynamic_cast<engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT> *>(other);
        if (o)
        {
            this->CorePointer = o->CorePointer;
            this->CoreTupleSize = o->CoreTupleSize;
            this->CoreValueSize = o->CoreValueSize;
            this->VoidPointerBuffer = o->VoidPointerBuffer;
            this->ExtraBuffer = o->ExtraBuffer;
            this->Size = o->Size;
            this->MaxId = o->MaxId;
            this->SetLookupTable(nullptr);
            this->Modified();
        }
    }
}

template<class CoreValueTypeT, typename BaseValueTypeT>
void engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>::InsertTuples(vtkIdType dstStart, vtkIdType n, vtkIdType srcStart, vtkAbstractArray* source)
{
    if (dstStart < this->CoreTupleSize)
    {
        vtkErrorMacro("Cannot insert at the core address");
        return;
    }
    if (source->GetNumberOfComponents() != this->GetNumberOfComponents())
    {
        vtkErrorMacro("Different number of components");
        return;
    }
    vtkDataArray* dataArray = this->SafeDownCast(source);
    if (!dataArray)
    {
        vtkErrorMacro("Not a data array");
        return;
    }
    vtkIdType neededSize = dstStart + n;
    vtkIdType neededExtraSize = neededSize - this->CoreTupleSize;
    if (neededExtraSize > static_cast<vtkIdType>(this->ExtraBuffer->size())) this->AllocateTuples(neededSize);

    for (vtkIdType i = dstStart; i < neededSize; i++)
    {
        this->SetTuple(i, dataArray->GetTuple(srcStart));
        srcStart++;
    }
}

template<class CoreValueTypeT, typename BaseValueTypeT>
void engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>::InsertTuples(vtkIdList* dstIds, vtkIdList* srcIds, vtkAbstractArray* source)
{
    if (dstIds->GetNumberOfIds() != srcIds->GetNumberOfIds())
    {
        vtkErrorMacro("Lists with different number of ids");
        return;
    }
    if (source->GetNumberOfComponents() != this->GetNumberOfComponents())
    {
        vtkErrorMacro("Different number of components");
        return;
    }
    vtkDataArray* dataArray = this->SafeDownCast(source);
    if (!dataArray)
    {
        vtkErrorMacro("Not a data array");
        return;
    }

    vtkIdType nIds = dstIds->GetNumberOfIds();
    for (vtkIdType i = 0; i < nIds; i++)
    {
        this->InsertTuple(dstIds->GetId(i), dataArray->GetTuple(srcIds->GetId(i)));
    }
}

template<class CoreValueTypeT, typename BaseValueTypeT>
void engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>::InsertTuplesStartingAt(vtkIdType dstStart, vtkIdList* srcIds, vtkAbstractArray* source)
{
    if (dstStart < this->CoreTupleSize)
    {
        vtkErrorMacro("Cannot insert at the core address");
        return;
    }
    if (source->GetNumberOfComponents() != this->GetNumberOfComponents())
    {
        vtkErrorMacro("Different number of components");
        return;
    }
    vtkDataArray* dataArray = this->SafeDownCast(source);
    if (!dataArray)
    {
        vtkErrorMacro("Not a data array");
        return;
    }

    vtkIdType n = srcIds->GetNumberOfIds();
    vtkIdType neededSize = dstStart + n;
    vtkIdType neededExtraSize = neededSize - this->CoreTupleSize;
    if (neededExtraSize > static_cast<vtkIdType>(this->ExtraBuffer->size())) this->AllocateTuples(neededSize);

    for (vtkIdType i = 0; i < n; i++)
    {
        this->SetTuple(dstStart, dataArray->GetTuple(srcIds->GetId(i)));
        dstStart++;
    }
}

template<class CoreValueTypeT, typename BaseValueTypeT>
void engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>::SetNumberOfTuples(vtkIdType nTuples)
{
    if (nTuples < this->CoreTupleSize)
    {
        vtkErrorMacro("Cannot resize to less than the core size");
        return;
    }

    vtkIdType neededExtraSize = nTuples - this->CoreTupleSize;
    if (neededExtraSize == 0)
    {
        this->ExtraBuffer->clear();
    }
    else
    {
        this->ExtraBuffer->resize(neededExtraSize);
        this->ExtraBuffer->shrink_to_fit();
    }
    this->Size = nTuples * this->CoreNumberOfComponents;
    this->MaxId = this->Size - 1;
}

template<class CoreValueTypeT, typename BaseValueTypeT>
void* engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>::GetVoidPointer(vtkIdType pos)
{
    if ((this->MaxId + 1) == (this->CoreTupleSize * this->CoreNumberOfComponents))
    {
        return const_cast<CoreValueTypeT*>(&this->CorePointer->operator[](pos));
    }
    else
    {
        if (!this->VoidPointerBuffer.get())
        {
            this->InitializeVoidPointerVariable();
        }
        return &this->VoidPointerBuffer->operator[](pos);
    }
}

template<class CoreValueTypeT, typename BaseValueTypeT>
void engysFoamCoreDataArray<CoreValueTypeT, BaseValueTypeT>::InitializeVoidPointerVariable()
{
    this->VoidPointerBuffer = std::shared_ptr<Foam::List<CoreValueTypeT>>(new Foam::List<CoreValueTypeT>(this->Size));

    Foam::label index = 0;
    for (Foam::label i = 0; i < this->CoreTupleSize; i++)
    {
        this->VoidPointerBuffer->operator[](index++) = this->CorePointer->operator[](i);
    }
    for (size_t i = 0; i < this->ExtraBuffer->size(); i++)
    {
        this->VoidPointerBuffer->operator[](index++) = this->ExtraBuffer->operator[](i);
    }
    this->CorePointer = this->VoidPointerBuffer.get();
    this->ExtraBuffer->clear();
    this->CoreTupleSize = this->VoidPointerBuffer->size();
    this->CoreValueSize = this->CoreTupleSize * this->CoreNumberOfComponents;
}


// ************************************************************************* //
