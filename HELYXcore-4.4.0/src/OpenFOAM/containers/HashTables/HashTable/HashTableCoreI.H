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
    (c) 2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * helper methods  * * * * * * * * * * * * * * //

template<class IteratorType, class TableType>
inline IteratorType Foam::HashTableCore::iterator_begin
(
    TableType& table
)
{
    return IteratorType(table.begin());
}


template<class IteratorType, class TableType>
inline IteratorType Foam::HashTableCore::iterator_begin
(
    const TableType& table
)
{
    return IteratorType(table.begin());
}


template<class IteratorType, class TableType>
inline IteratorType Foam::HashTableCore::iterator_cbegin
(
    const TableType& table
)
{
    return IteratorType(table.cbegin());
}


// * * * * * * * * * * * * * const iterator pair * * * * * * * * * * * * * * //

template<class IteratorType, class TableType>
inline Foam::HashTableCore::const_iterator_pair<IteratorType, TableType>
::const_iterator_pair
(
    const TableType& tbl
)
:
    size_(tbl.size()),
    iter_(tbl.begin()),
    end_(tbl.end())
{}


template<class IteratorType, class TableType>
inline Foam::label
Foam::HashTableCore::const_iterator_pair<IteratorType, TableType>::size() const
{
    return size_;
}


template<class IteratorType, class TableType>
inline bool
Foam::HashTableCore::const_iterator_pair<IteratorType, TableType>::empty() const
{
    return !size_;
}


template<class IteratorType, class TableType>
inline IteratorType Foam::HashTableCore::const_iterator_pair
<
    IteratorType,
    TableType
>::begin() const
{
    return iter_;
}


template<class IteratorType, class TableType>
inline IteratorType Foam::HashTableCore::const_iterator_pair
<
    IteratorType,
    TableType
>::cbegin() const
{
    return iter_;
}


template<class IteratorType, class TableType>
inline const IteratorType Foam::HashTableCore::const_iterator_pair
<
    IteratorType,
    TableType
>::end() const
{
    return end_;
}


template<class IteratorType, class TableType>
inline const IteratorType Foam::HashTableCore::const_iterator_pair
<
    IteratorType,
    TableType
>::cend() const
{
    return end_;
}


// ************************************************************************* //
