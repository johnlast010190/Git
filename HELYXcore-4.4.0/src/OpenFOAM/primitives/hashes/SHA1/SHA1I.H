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

#include "primitives/hashes/SHA1/SHA1.H"
#include <cstring>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::SHA1::SHA1()
{
    clear();
}


inline Foam::SHA1::SHA1(const std::string& str)
{
    clear();
    append(str);
}


inline Foam::SHA1::SHA1(const char* str)
{
    clear();
    append(str);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::SHA1& Foam::SHA1::append(const char* data, size_t len)
{
    processBytes(data, len);
    return *this;
}


inline Foam::SHA1& Foam::SHA1::append(const std::string& str)
{
    processBytes(str.data(), str.size());
    return *this;
}


inline Foam::SHA1& Foam::SHA1::append(const char* str)
{
    if (str)
    {
        processBytes(str, strlen(str));
    }
    return *this;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::SHA1::operator==(const SHA1& rhs) const
{
    return this->digest() == rhs.digest();
}


inline bool Foam::SHA1::operator==(const SHA1Digest& rhs) const
{
    return this->digest() == rhs;
}


inline bool Foam::SHA1::operator==(const std::string& hexdigits) const
{
    return this->digest() == hexdigits;
}


inline bool Foam::SHA1::operator==(const char* hexdigits) const
{
    return this->digest() == hexdigits;
}


inline bool Foam::SHA1::operator!=(const SHA1& rhs) const
{
    return !this->operator==(rhs);
}


inline bool Foam::SHA1::operator!=(const SHA1Digest& rhs) const
{
    return !this->operator==(rhs);
}


inline bool Foam::SHA1::operator!=(const std::string& rhs) const
{
    return !this->operator==(rhs);
}


inline bool Foam::SHA1::operator!=(const char* rhs) const
{
    return !this->operator==(rhs);
}


inline Foam::SHA1::operator Foam::SHA1Digest() const
{
    return digest();
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<(Ostream& os, const SHA1& sha)
{
    return os  << sha.digest();
}


// ************************************************************************* //
