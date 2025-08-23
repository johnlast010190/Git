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

\*---------------------------------------------------------------------------*/

#include "containers/Lists/ListOps/ListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::SortableList<T>::SortableList()
{}


template<class T>
Foam::SortableList<T>::SortableList(const UList<T>& values)
:
    List<T>(values)
{
    sort();
}


template<class T>
Foam::SortableList<T>::SortableList(List<T>&& values)
:
    List<T>(std::move(values))
{
    sort();
}


template<class T>
inline Foam::SortableList<T>::SortableList(const label size)
:
    List<T>(size)
{}


template<class T>
inline Foam::SortableList<T>::SortableList(const label size, const T& val)
:
    List<T>(size, val)
{}


template<class T>
Foam::SortableList<T>::SortableList(const SortableList<T>& lst)
:
    List<T>(lst),
    indices_(lst.indices())
{}


template<class T>
template<class InputIterator>
inline Foam::SortableList<T>::SortableList
(
    InputIterator begIter,
    InputIterator endIter
)
:
    List<T>(begIter, endIter)
{
    sort();
}


template<class T>
Foam::SortableList<T>::SortableList(SortableList<T>&& lst)
:
    List<T>(std::move(lst)),
    indices_(std::move(lst.indices()))
{}


template<class T>
Foam::SortableList<T>::SortableList(std::initializer_list<T> values)
:
    List<T>(values)
{
    sort();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class T>
void Foam::SortableList<T>::clear()
{
    List<T>::clear();
    indices_.clear();
}


template<class T>
Foam::List<T>& Foam::SortableList<T>::shrink()
{
    indices_.clear();
    return static_cast<List<T>&>(*this);
}


template<class T>
void Foam::SortableList<T>::sort()
{
    sortedOrder(*this, indices_);

    List<T> lst(this->size());
    forAll(indices_, i)
    {
        lst[i] = this->operator[](indices_[i]);
    }

    List<T>::transfer(lst);
}


template<class T>
void Foam::SortableList<T>::reverseSort()
{
    sortedOrder(*this, indices_, typename UList<T>::greater(*this));

    List<T> lst(this->size());
    forAll(indices_, i)
    {
        lst[i] = this->operator[](indices_[i]);
    }

    List<T>::transfer(lst);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline void Foam::SortableList<T>::operator=(const T& val)
{
    UList<T>::operator=(val);
}


template<class T>
inline void Foam::SortableList<T>::operator=(const UList<T>& lst)
{
    List<T>::operator=(lst);
    indices_.clear();
}


template<class T>
inline void Foam::SortableList<T>::operator=(const SortableList<T>& lst)
{
    List<T>::operator=(lst);
    indices_ = lst.indices();
}


template<class T>
inline void Foam::SortableList<T>::operator=(SortableList<T>&& lst)
{
    List<T>::operator=(std::move(lst));
    indices_ = std::move(lst.indices());
}



template<class T>
inline void Foam::SortableList<T>::operator=(std::initializer_list<T> lst)
{
    List<T>::operator=(lst);
    sort();
}


// ************************************************************************* //
