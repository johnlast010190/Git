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
    (c) 2011-2014 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "db/IOstreams/Sstreams/prefixOSstream.H"
#include "db/IOstreams/token/token.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

inline void Foam::prefixOSstream::checkWritePrefix()
{
    if (printPrefix_ && prefix_.size())
    {
        OSstream::write(prefix_.c_str());
        printPrefix_ = false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::prefixOSstream::prefixOSstream
(
    ostream& os,
    const string& name,
    streamFormat format,
    versionNumber version,
    compressionType compression
)
:
    OSstream(os, name, format, version, compression),
    printPrefix_(true),
    prefix_("")
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::prefixOSstream::print(Ostream& os) const
{
    os  << "prefixOSstream ";
    OSstream::print(os);
}


Foam::Ostream& Foam::prefixOSstream::write(const token& t)
{
    if (t.type() == token::VERBATIMSTRING)
    {
        write(char(token::HASH));
        write(char(token::BEGIN_BLOCK));
        writeQuoted(t.stringToken(), false);
        write(char(token::HASH));
        write(char(token::END_BLOCK));
    }
    else if (t.type() == token::VARIABLE)
    {
        writeQuoted(t.stringToken(), false);
    }
    return *this;
}


Foam::Ostream& Foam::prefixOSstream::write(const char c)
{
    checkWritePrefix();
    OSstream::write(c);

    if (c == token::NL)
    {
        printPrefix_ = true;
    }

    return *this;
}


Foam::Ostream& Foam::prefixOSstream::write(const char* str)
{
    checkWritePrefix();
    OSstream::write(str);

    size_t len = strlen(str);
    if (len && str[len-1] == token::NL)
    {
        printPrefix_ = true;
    }

    return *this;
}


Foam::Ostream& Foam::prefixOSstream::write(const word& val)
{
    checkWritePrefix();
    return OSstream::write(val);
}


Foam::Ostream& Foam::prefixOSstream::write(const string& val)
{
    checkWritePrefix();
    return OSstream::write(val);
}


Foam::Ostream& Foam::prefixOSstream::writeQuoted
(
    const std::string& val,
    const bool quoted
)
{
    checkWritePrefix();
    return OSstream::writeQuoted(val, quoted);
}


Foam::Ostream& Foam::prefixOSstream::write(const int32_t val)
{
    checkWritePrefix();
    return OSstream::write(val);
}


Foam::Ostream& Foam::prefixOSstream::write(const int64_t val)
{
    checkWritePrefix();
    return OSstream::write(val);
}


Foam::Ostream& Foam::prefixOSstream::write(const floatScalar val)
{
    checkWritePrefix();
    return OSstream::write(val);
}


Foam::Ostream& Foam::prefixOSstream::write(const doubleScalar val)
{
    checkWritePrefix();
    return OSstream::write(val);
}


Foam::Ostream& Foam::prefixOSstream::write
(
    const char* buf,
    std::streamsize count
)
{
    checkWritePrefix();
    return OSstream::write(buf, count);
}


void Foam::prefixOSstream::indent()
{
    checkWritePrefix();
    OSstream::indent();
}

// ************************************************************************* //
