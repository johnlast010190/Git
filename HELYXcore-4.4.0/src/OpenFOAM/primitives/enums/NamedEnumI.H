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
    (c) 2011-2018 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EnumType, int nEnum>
inline Foam::label Foam::NamedEnum<EnumType, nEnum>::size() const
{
    return lookup_.size();
}


template<class EnumType, int nEnum>
inline Foam::wordList Foam::NamedEnum<EnumType, nEnum>::toc() const
{
    return lookup_.toc();
}


template<class EnumType, int nEnum>
inline Foam::wordList Foam::NamedEnum<EnumType, nEnum>::sortedToc() const
{
    return lookup_.sortedToc();
}


template<class EnumType, int nEnum>
inline bool Foam::NamedEnum<EnumType, nEnum>::found
(
    const word& enumName
) const
{
    return lookup_.found(enumName);
}


template<class EnumType, int nEnum>
inline bool Foam::NamedEnum<EnumType, nEnum>::hasEnum
(
    const word& enumName
) const
{
    return lookup_.found(enumName);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EnumType, int nEnum>
inline const EnumType Foam::NamedEnum<EnumType, nEnum>::operator[]
(
    const word& name
) const
{
    return EnumType(lookup_[name]);
}


template<class EnumType, int nEnum>
inline const char* Foam::NamedEnum<EnumType, nEnum>::operator[]
(
    const EnumType e
) const
{
    unsigned int ue = unsigned(e);

    if (ue < nEnum)
    {
        return names[ue];
    }
    else
    {
        FatalErrorInFunction
            << "names array index " << ue << " out of range 0-"
            << nEnum - 1
            << exit(FatalError);

        return names[0];
    }
}


// ************************************************************************* //
