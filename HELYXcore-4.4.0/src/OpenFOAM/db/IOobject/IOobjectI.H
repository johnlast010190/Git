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
    (c) 2011-2020 OpenFOAM Foundation
    (c) 2017-2019 OpenCFD Ltd
    (c) 2010-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class StringType>
inline Foam::word Foam::IOobject::groupName(StringType name, const word& group)
{
    if (group.empty())
    {
        return name;
    }
    else
    {
        return name + ('.' + group);
    }
}


template<class Name>
inline Foam::word Foam::IOobject::modelName(Name name, const word& model)
{
    if (model != word::null)
    {
        return name + ('-' + model);
    }
    else
    {
        return name;
    }
}


inline Foam::word Foam::IOobject::modelName(const char* name) const
{
    return type() + '-' + name;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word& Foam::IOobject::name() const
{
    return name_;
}


inline const Foam::word& Foam::IOobject::headerClassName() const
{
    return headerClassName_;
}


inline Foam::word& Foam::IOobject::headerClassName()
{
    return headerClassName_;
}


inline const Foam::string& Foam::IOobject::note() const
{
    return note_;
}


inline Foam::string& Foam::IOobject::note()
{
    return note_;
}


inline bool Foam::IOobject::registerObject() const
{
    return registerObject_;
}


inline bool& Foam::IOobject::registerObject()
{
    return registerObject_;
}


inline Foam::IOobject::readOption Foam::IOobject::readOpt() const
{
    return rOpt_;
}


inline Foam::IOobject::readOption& Foam::IOobject::readOpt()
{
    return rOpt_;
}


inline Foam::IOobject::writeOption Foam::IOobject::writeOpt() const
{
    return wOpt_;
}


inline Foam::IOobject::writeOption& Foam::IOobject::writeOpt()
{
    return wOpt_;
}


inline const Foam::fileName& Foam::IOobject::instance() const
{
    return instance_;
}


inline Foam::fileName& Foam::IOobject::instance()
{
    return instance_;
}


inline const Foam::fileName& Foam::IOobject::local() const
{
    return local_;
}


inline Foam::fileName Foam::IOobject::objectPath() const
{
#if defined(WIN64) || defined(WIN32)
    return path()/uniqueFileName();
#else
    return path()/name();
#endif
}


inline bool Foam::IOobject::good() const
{
    return objState_ == GOOD;
}


inline bool Foam::IOobject::bad() const
{
    return objState_ == BAD;
}


inline bool Foam::IOobject::isHeaderClassName(const word& clsName) const
{
    return (clsName == headerClassName_);
}


template<class Type>
inline bool Foam::IOobject::isHeaderClassName() const
{
    return (Type::typeName == headerClassName_);
}

// ************************************************************************* //
