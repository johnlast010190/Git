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

#include "db/IOobjects/IOMap/IOMap.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::IOMap<T>::IOMap(const IOobject& io)
:
    regIOobject(io)
{
    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        // For if MUST_READ_IF_MODIFIED
        addWatch();

        readStream(typeName) >> *this;
        close();
    }
}

template<class T>
Foam::IOMap<T>::IOMap(const IOobject& io, const label size)
:
    regIOobject(io)
{
    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        // For if MUST_READ_IF_MODIFIED
        addWatch();

        readStream(typeName) >> *this;
        close();
    }
    else
    {
        Map<T>::setSize(size);
    }
}


template<class T>
Foam::IOMap<T>::IOMap(const IOobject& io, const Map<T>& map)
:
    regIOobject(io)
{
    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        // For if MUST_READ_IF_MODIFIED
        addWatch();

        readStream(typeName) >> *this;
        close();
    }
    else
    {
        Map<T>::operator=(map);
    }
}


template<class T>
Foam::IOMap<T>::IOMap(const IOobject& io, Map<T>&& map)
:
    regIOobject(io),
    Map<T>(std::move(map))
{
    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        // For if MUST_READ_IF_MODIFIED
        addWatch();

        readStream(typeName) >> *this;
        close();
    }
}


template<class T>
Foam::IOMap<T>::IOMap(IOMap<T>&& map)
:
    regIOobject(std::move(map)),
    Map<T>(std::move(map))
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class T>
Foam::IOMap<T>::~IOMap()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
bool Foam::IOMap<T>::writeData(Ostream& os) const
{
    return (os << *this).good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
void Foam::IOMap<T>::operator=(const IOMap<T>& rhs)
{
    Map<T>::operator=(rhs);
}


template<class T>
void Foam::IOMap<T>::operator=(IOMap<T>&& rhs)
{
    Map<T>::operator=(std::move(rhs));
}


template<class T>
void Foam::IOMap<T>::operator=(const Map<T>& rhs)
{
    Map<T>::operator=(rhs);
}


template<class T>
void Foam::IOMap<T>::operator=(Map<T>&& rhs)
{
    Map<T>::operator=(std::move(rhs));
}


// ************************************************************************* //
