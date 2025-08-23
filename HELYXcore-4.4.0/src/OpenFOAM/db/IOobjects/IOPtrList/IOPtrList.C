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

#include "db/IOobjects/IOPtrList/IOPtrList.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
template<class INew>
Foam::IOPtrList<T>::IOPtrList(const IOobject& io, const INew& inewt)
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

        PtrList<T>::read(readStream(typeName), inewt);
        close();
    }
}


template<class T>
Foam::IOPtrList<T>::IOPtrList(const IOobject& io)
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

        PtrList<T>::read(readStream(typeName), INew<T>());
        close();
    }
}


template<class T>
Foam::IOPtrList<T>::IOPtrList(const IOobject& io, const label s)
:
    regIOobject(io),
    PtrList<T>(s)
{
    if (io.readOpt() != IOobject::NO_READ)
    {
        FatalErrorInFunction
            << "NO_READ must be set if specifying size" << nl
            << exit(FatalError);
    }
}


template<class T>
Foam::IOPtrList<T>::IOPtrList(const IOobject& io, const PtrList<T>& list)
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

        PtrList<T>::read(readStream(typeName), INew<T>());
        close();
    }
    else
    {
        PtrList<T>::operator=(list);
    }
}


template<class T>
Foam::IOPtrList<T>::IOPtrList(const IOobject& io, PtrList<T>&& list)
:
    regIOobject(io),
    PtrList<T>(std::move(list))
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

        PtrList<T>::read(readStream(typeName), INew<T>());
        close();
    }
}


template<class T>
Foam::IOPtrList<T>::IOPtrList(IOPtrList<T>&& list)
:
    regIOobject(std::move(list)),
    PtrList<T>(std::move(list))
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class T>
Foam::IOPtrList<T>::~IOPtrList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
bool Foam::IOPtrList<T>::writeData(Ostream& os) const
{
    return (os << *this).good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
void Foam::IOPtrList<T>::operator=(const IOPtrList<T>& rhs)
{
    PtrList<T>::operator=(rhs);
}


template<class T>
void Foam::IOPtrList<T>::operator=(IOPtrList<T>&& rhs)
{
    PtrList<T>::operator=(std::move(rhs));
}


// ************************************************************************* //
