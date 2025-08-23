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

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function1/Table/Table.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::Table<Type>::Table
(
    const word& entryName,
    const dictionary& dict
)
:
    TableBase<Type>(entryName, dict)
{
    Istream& is(dict.lookup(entryName));
    word entryType(is);
    is  >> this->table_;
    TableBase<Type>::check();
}


template<class Type>
Foam::Function1Types::Table<Type>::Table(const Table<Type>& tbl)
:
    TableBase<Type>(tbl)
{}


template<class Type>
Foam::Function1Types::Table<Type>::Table
(
    const word& entryName,
    const List<Tuple2<scalar, Type>>& tbl
)
:
    TableBase<Type>(entryName, dictionary())
{
    this->table_ = tbl;
    TableBase<Type>::check();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::Table<Type>::~Table()
{}


// ************************************************************************* //
