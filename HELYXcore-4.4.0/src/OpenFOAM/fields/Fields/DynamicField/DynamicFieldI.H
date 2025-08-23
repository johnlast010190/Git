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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::DynamicField()
:
    Field<T>(0),
    capacity_(Field<T>::size())
{}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::DynamicField
(
    const label nElem
)
:
    Field<T>(nElem),
    capacity_(Field<T>::size())
{
    // we could also enforce SizeInc granularity when (!SizeMult || !SizeDiv)
    Field<T>::size(0);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::DynamicField
(
    const UList<T>& lst
)
:
    Field<T>(lst),
    capacity_(Field<T>::size())
{}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::DynamicField
(
    List<T>&& lst
)
:
    Field<T>(std::move(lst)),
    capacity_(Field<T>::size())
{}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::DynamicField
(
    Field<T>&& lst
)
:
    Field<T>(std::move(lst)),
    capacity_(Field<T>::size())
{}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::DynamicField
(
    const UList<T>& mapF,
    const labelList& mapAddressing
)
:
    Field<T>(mapF, mapAddressing),
    capacity_(Field<T>::size())
{}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::DynamicField
(
    const UList<T>& mapF,
    const labelListList& mapAddressing,
    const scalarListList& weights
)
:
    Field<T>(mapF, mapAddressing, weights),
    capacity_(Field<T>::size())
{}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::DynamicField
(
    const DynamicField<T, SizeInc, SizeMult, SizeDiv>& lst
)
:
    Field<T>(std::move(lst)),
    capacity_(lst.capacity_)
{}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::DynamicField
(
    DynamicField<T, SizeInc, SizeMult, SizeDiv>&& lst
)
:
    Field<T>(std::move(lst)),
    capacity_(lst.capacity_)
{
    lst.capacity_ = 0;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::label Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::capacity()
const
{
    return capacity_;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::setCapacity
(
    const label nElem
)
{
    label nextFree = Field<T>::size();
    capacity_ = nElem;

    if (nextFree > capacity_)
    {
        // truncate addressed sizes too
        nextFree = capacity_;
    }
    // we could also enforce SizeInc granularity when (!SizeMult || !SizeDiv)

    Field<T>::setSize(capacity_);
    Field<T>::size(nextFree);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::reserve
(
    const label nElem
)
{
    // allocate more capacity?
    if (nElem > capacity_)
    {
        capacity_ =
            max
            (
                nElem,
                label(SizeInc + capacity_ * SizeMult / SizeDiv)
            );

        // adjust allocated size, leave addressed size untouched
        label nextFree = Field<T>::size();
        Field<T>::setSize(capacity_);
        Field<T>::size(nextFree);
    }
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::setSize
(
    const label nElem
)
{
    // allocate more capacity?
    if (nElem > capacity_)
    {
        capacity_ =
            max
            (
                nElem,
                label(SizeInc + capacity_ * SizeMult / SizeDiv)
            );

        Field<T>::setSize(capacity_);
    }

    // adjust addressed size
    Field<T>::size(nElem);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::setSize
(
    const label nElem,
    const T& t
)
{
    label nextFree = Field<T>::size();
    setSize(nElem);

    // set new elements to constant value
    while (nextFree < nElem)
    {
        this->operator[](nextFree++) = t;
    }
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::resize
(
    const label nElem
)
{
    this->setSize(nElem);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::resize
(
    const label nElem,
    const T& t
)
{
    this->setSize(nElem, t);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::clear()
{
    Field<T>::size(0);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::clearStorage()
{
    Field<T>::clear();
    capacity_ = 0;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>&
Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::shrink()
{
    label nextFree = Field<T>::size();
    if (capacity_ > nextFree)
    {
        // use the full list when resizing
        Field<T>::size(capacity_);

        // the new size
        capacity_ = nextFree;
        Field<T>::setSize(capacity_);
        Field<T>::size(nextFree);
    }
    return *this;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>&
Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::append
(
    const T& t
)
{
    const label elemI = List<T>::size();
    setSize(elemI + 1);

    this->operator[](elemI) = t;
    return *this;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>&
Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::append
(
    const UList<T>& lst
)
{
    if (this == &lst)
    {
        FatalErrorInFunction
            << "attempted appending to self" << abort(FatalError);
    }

    label nextFree = List<T>::size();
    setSize(nextFree + lst.size());

    forAll(lst, elemI)
    {
        this->operator[](nextFree++) = lst[elemI];
    }
    return *this;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline T Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::remove()
{
    const label elemI = List<T>::size() - 1;

    if (elemI < 0)
    {
        FatalErrorInFunction
            << "List is empty" << abort(FatalError);
    }

    const T& val = List<T>::operator[](elemI);

    List<T>::size(elemI);

    return val;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline T& Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::operator()
(
    const label elemI
)
{
    if (elemI >= Field<T>::size())
    {
        setSize(elemI + 1);
    }

    return this->operator[](elemI);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::operator=
(
    const T& t
)
{
    UList<T>::operator=(t);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::operator=
(
    const DynamicField<T, SizeInc, SizeMult, SizeDiv>& lst
)
{
    if (this == &lst)
    {
        FatalErrorInFunction
            << "attempted assignment to self" << abort(FatalError);
    }

    if (capacity_ >= lst.size())
    {
        // can copy w/o reallocating, match initial size to avoid reallocation
        Field<T>::size(lst.size());
        Field<T>::operator=(lst);
    }
    else
    {
        // make everything available for the copy operation
        Field<T>::size(capacity_);

        Field<T>::operator=(lst);
        capacity_ = Field<T>::size();
    }
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::operator=
(
    DynamicField<T, SizeInc, SizeMult, SizeDiv>&& lst
)
{
    if (this == &lst)
    {
        FatalErrorInFunction
            << "attempted assignment to self" << abort(FatalError);
    }

    Field<T>::operator=(std::move(lst));
    capacity_ = lst.capacity_;
    lst.capacity_ = 0;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::operator=
(
    const UList<T>& lst
)
{
    if (capacity_ >= lst.size())
    {
        // can copy w/o reallocating, match initial size to avoid reallocation
        Field<T>::size(lst.size());
        Field<T>::operator=(lst);
    }
    else
    {
        // make everything available for the copy operation
        Field<T>::size(capacity_);

        Field<T>::operator=(lst);
        capacity_ = Field<T>::size();
    }
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::operator=
(
    List<T>&& lst
)
{
    if (this == &lst)
    {
        FatalErrorInFunction
            << "attempted assignment to self" << abort(FatalError);
    }

    Field<T>::operator=(std::move(lst));
    capacity_ = Field<T>::size();
}


// ************************************************************************* //
