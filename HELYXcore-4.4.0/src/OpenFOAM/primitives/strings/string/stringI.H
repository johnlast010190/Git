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
    (c) 2018-2019 OpenCFD Ltd.
    (c) 2011-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

inline std::string::size_type Foam::string::find_ext() const
{
    const size_type i = find_last_of("./");

    if (i == npos || i == 0 || operator[](i) == '/')
    {
        return npos;
    }
    else
    {
        return i;
    }
}


inline bool Foam::string::hasExt() const
{
    return (find_ext() != npos);
}


inline bool Foam::string::removePath()
{
    const auto i = rfind('/');

    if (npos == i)
    {
        return false;
    }

    this->erase(0, i+1);
    return true;
}


inline bool Foam::string::removeExt()
{
    const size_type i = find_ext();

    if (i == npos)
    {
        return false;
    }
    else
    {
        this->resize(i);
        return true;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::string::string()
{}


inline Foam::string::string(const std::string& str)
:
    std::string(str)
{}


inline Foam::string::string(const string& str)
:
    std::string(str)
{}


inline Foam::string::string(string&& str)
:
    std::string(std::move(str))
{}


inline Foam::string::string(const char* str)
:
    std::string(str)
{}


inline Foam::string::string(const char* str, const size_type len)
:
    std::string(str, len)
{}


inline Foam::string::string(const char c)
:
    std::string(1, c)
{}


inline Foam::string::string(const size_type len, const char c)
:
    std::string(len, c)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class String>
inline bool Foam::string::valid(const std::string& str)
{
    for (auto iter = str.cbegin(); iter != str.cend(); ++iter)
    {
        if (!String::valid(*iter))
        {
            return false;
        }
    }

    return true;
}


template<class String>
inline bool Foam::string::stripInvalid(std::string& str)
{
    if (!valid<String>(str))
    {
        size_type nChar = 0;
        iterator outIter = str.begin();

        for (auto iter = str.cbegin(); iter != str.cend(); ++iter)
        {
            const char c = *iter;

            if (String::valid(c))
            {
                *outIter = c;
                ++outIter;
                ++nChar;
            }
        }

        str.resize(nChar);

        return true;
    }

    return false;
}


template<class String>
inline bool Foam::string::meta(const std::string& str, const char quote)
{
    int escaped = 0;
    for (auto iter = str.cbegin(); iter != str.cend(); ++iter)
    {
        const char c = *iter;
        if (quote && c == quote)
        {
            escaped ^= 1;  // toggle state
        }
        else if (escaped)
        {
            escaped = 0;
        }
        else if (String::meta(c))
        {
            return true;
        }
    }
    return false;
}


template<class String>
inline Foam::string
Foam::string::quotemeta(const std::string& str, const char quote)
{
    if (!quote)
    {
        return str;
    }

    string sQuoted;
    sQuoted.reserve(2*str.size());

    int escaped = 0;
    for (auto iter = str.cbegin(); iter != str.cend(); ++iter)
    {
        const char c = *iter;
        if (c == quote)
        {
            escaped ^= 1;  // toggle state
        }
        else if (escaped)
        {
            escaped = 0;
        }
        else if (String::meta(c))
        {
            sQuoted += quote;
        }

        sQuoted += c;
    }

    sQuoted.resize(sQuoted.size());

    return sQuoted;
}


template<class String>
inline String Foam::string::validate(const std::string& str)
{
    String out;
    out.resize(str.size());

    size_type count = 0;
    for (auto iter = str.cbegin(); iter != str.cend(); ++iter)
    {
        const char c = *iter;
        if (String::valid(c))
        {
            out[count++] = c;
        }
    }

    out.resize(count);

    return out;
}


inline bool Foam::string::match(const std::string& text) const
{
    return !compare(text);  // Always compare as literal string
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::string::operator()(const std::string& text) const
{
    return !compare(text);  // Always compare as literal string
}


inline Foam::string Foam::string::operator()
(
    const size_type i,
    const size_type n
) const
{
    return substr(i, n);
}


inline Foam::string Foam::string::operator()(const size_type n) const
{
    return substr(0, n);
}


inline unsigned Foam::string::hash::operator()
(
    const string& str,
    unsigned seed
) const
{
    return Hasher(str.data(), str.size(), seed);
}


void Foam::string::operator=(const string& str)
{
    std::string::operator=(str);
}


void Foam::string::operator=(string&& str)
{
    std::string::operator=(std::move(str));
}


// ************************************************************************* //
