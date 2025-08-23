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
    (c) Creative Fields, Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Description
    An IOLongList of a given type is a list which supports automated
    input and output.

\*---------------------------------------------------------------------------*/

#include "utilities/containers/IOLongList/IOLongList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T, label Offset>
IOLongList<T, Offset>::IOLongList(const IOobject& io)
:
    regIOobject(io),
    LongList<T, Offset>()
{
    if
    (
        io.readOpt() == IOobject::MUST_READ
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }
}


template<class T, label Offset>
IOLongList<T, Offset>::IOLongList
(
    const IOobject& io,
    const label size
)
:
    regIOobject(io),
    LongList<T, Offset>(size)
{}


template<class T, label Offset>
IOLongList<T, Offset>::IOLongList
(
    const IOobject& io,
    const LongList<T, Offset>& list
)
:
    regIOobject(io),
    LongList<T, Offset>()
{
    if (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    {
        readStream(typeName) >> *this;
        close();
    }

    LongList<T, Offset>::operator=(list);
}


template<class T, label Offset>
void IOLongList<T, Offset>::operator=
(
    const IOLongList<T, Offset>& rhs
)
{
    LongList<T, Offset>::operator=(rhs);
}


template<class T, label Offset>
void IOLongList<T, Offset>::operator=
(
    const LongList<T, Offset>& rhs
)
{
    LongList<T, Offset>::operator=(rhs);
}


template<class T, label Offset>
bool IOLongList<T, Offset>::writeData(Ostream& os) const
{
    return (os << *this).good();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
