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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


inline Foam::IndirectListAddressing::IndirectListAddressing
(
    const labelUList& addr
)
:
    addressing_(addr)
{}


inline Foam::IndirectListAddressing::IndirectListAddressing
(
    List<label>&& addr
)
:
    addressing_(std::move(addr))
{}


template<class T>
inline Foam::IndirectList<T>::IndirectList
(
    const UList<T>& completeList,
    const labelUList& addr
)
:
    IndirectListAddressing(addr),
    UIndirectList<T>
    (
        completeList,
        IndirectListAddressing::addressing()
    )
{}


template<class T>
inline Foam::IndirectList<T>::IndirectList
(
    const UList<T>& completeList,
    List<label>&& addr
)
:
    IndirectListAddressing(std::move(addr)),
    UIndirectList<T>
    (
        completeList,
        IndirectListAddressing::addressing()
    )
{}


template<class T>
inline Foam::IndirectList<T>::IndirectList
(
    const IndirectList<T>& lst
)
:
    IndirectListAddressing(lst.addressing()),
    UIndirectList<T>
    (
        lst.completeList(),
        IndirectListAddressing::addressing()
    )
{}


template<class T>
inline Foam::IndirectList<T>::IndirectList
(
    const UIndirectList<T>& lst
)
:
    IndirectListAddressing(lst.addressing()),
    UIndirectList<T>
    (
        lst.completeList(),
        IndirectListAddressing::addressing()
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::List<Foam::label>&
Foam::IndirectListAddressing::addressing() const
{
    return addressing_;
}


inline void Foam::IndirectListAddressing::resetAddressing
(
    const labelUList& addr
)
{
    addressing_ = addr;
}


template<class T>
inline void Foam::IndirectList<T>::shallowCopy(const IndirectList<T>& a)
{
    UIndirectList<T>::operator=(a);
}


// ************************************************************************* //
