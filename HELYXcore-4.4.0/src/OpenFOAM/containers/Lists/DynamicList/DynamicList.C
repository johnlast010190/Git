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
    (c) 2011 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "containers/Lists/DynamicList/DynamicList.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::DynamicList(Istream& is)
:
    List<T>(is),
    capacity_(List<T>::size())
{}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const DynamicList<T, SizeInc, SizeMult, SizeDiv>& lst
)
{
    os << static_cast<const List<T>&>(lst);
    return os;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Foam::Istream& Foam::operator>>
(
    Istream& is,
    DynamicList<T, SizeInc, SizeMult, SizeDiv>& lst
)
{
    is >> static_cast<List<T>&>(lst);
    lst.capacity_ = lst.List<T>::size();

    return is;
}


// ************************************************************************* //
