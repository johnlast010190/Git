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
    (c) 2016-2019 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include <iostream>  // For std::cerr

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::fileName::stripInvalid()
{
    // skip stripping unless debug is active to avoid
    // costly operations
    if (debug && string::stripInvalid<fileName>(*this))
    {
        std::cerr
            << "fileName::stripInvalid() called for invalid fileName "
            << this->c_str() << std::endl;

        if (debug > 1)
        {
            std::cerr
                << "    For debug level (= " << debug
                << ") > 1 this is considered fatal" << std::endl;
            std::abort();
        }

        removeRepeated('/');
        removeTrailing('/');
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::fileName::fileName()
:
    string()
{}


inline Foam::fileName::fileName(const fileName& fn)
:
    string(fn)
{}


inline Foam::fileName::fileName(fileName&& fn)
:
    string(std::move(fn))
{}


inline Foam::fileName::fileName(const word& w)
:
    string(w)
#if defined(WIN64) || defined(WIN32)
{
    toUnixPath();
    stripInvalid();
}
#else
{}
#endif

inline Foam::fileName::fileName(const string& s, const bool doStripInvalid)
:
    string(s)
{
    if (doStripInvalid)
    {
#if defined(WIN64) || defined(WIN32)
        toUnixPath();
        #endif
        stripInvalid();
    }
}


inline Foam::fileName::fileName(const std::string& s, const bool doStripInvalid)
:
    string(s)
{
    if (doStripInvalid)
    {
#if defined(WIN64) || defined(WIN32)
        toUnixPath();
        #endif
        stripInvalid();
    }
}


inline Foam::fileName::fileName(const char* s, const bool doStripInvalid)
:
    string(s)
{
    if (doStripInvalid)
    {
#if defined(WIN64) || defined(WIN32)
        toUnixPath();
        #endif
        stripInvalid();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::fileName::valid(char c)
{
    return
    (
        !isspace(c)
     && c != '"'   // string quote
     && c != '\''  // string quote
    );
}

inline bool Foam::fileName::isAbsolute(const std::string& str)
{
    return
    (
        (!str.empty() && str.front() == '/')  // ie, str.starts_with('/')
        #if defined(WIN32) || defined(WIN64)
        ||
        (
            // Eg, d:/path or \\machine/path
            (str.length() > 1)
         && (str[1] == ':' || (str[0] == '\\' && str[1] == '\\'))
        )
        #endif
    );
}


inline bool Foam::fileName::isAbsolute() const
{
    return isAbsolute(*this);
}


inline bool Foam::fileName::hasExt() const
{
    return string::hasExt();
}

inline bool Foam::fileName::removePath()
{
    return string::removePath();
}


inline bool Foam::fileName::removeExt()
{
    return string::removeExt();
}


// ************************************************************************* //
