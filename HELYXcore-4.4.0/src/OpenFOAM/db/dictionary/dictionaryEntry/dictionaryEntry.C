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
    (c) 2010-2012 Engys Ltd.
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2018 OpenCFD Ltd.
\*---------------------------------------------------------------------------*/

#include "db/dictionary/dictionaryEntry/dictionaryEntry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dictionaryEntry::dictionaryEntry
(
    const keyType& key,
    const dictionary& parentDict,
    const dictionary& dict
)
:
    entry(key),
    dictionary(parentDict, dict)
{}


Foam::dictionaryEntry::dictionaryEntry
(
    const dictionary& parentDict,
    const dictionaryEntry& dictEnt
)
:
    entry(dictEnt),
    dictionary(parentDict, dictEnt)
{}

Foam::dictionaryEntry::dictionaryEntry
(
)
:
    entry(word("null")),
    dictionary(dictionary::null, dictionary())
{}

/*
Foam::dictionaryEntry::dictionaryEntry
(
    const dictionaryEntry& dictEntry
)
:
    refCount(),
    entry(dictEntry),
    dictionary(dictEntry.dict().parent(), dictEntry.dict())
{}*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::dictionaryEntry::startLineNumber() const
{
    if (size())
    {
        return first()->startLineNumber();
    }
    else
    {
        return -1;
    }
}


Foam::label Foam::dictionaryEntry::endLineNumber() const
{
    if (size())
    {
        return last()->endLineNumber();
    }
    else
    {
        return -1;
    }
}


Foam::ITstream& Foam::dictionaryEntry::stream() const
{
    FatalIOErrorInFunction(*this)
        << "Attempt to return dictionary entry as a primitive"
        << abort(FatalIOError);

    return lookup("");
}


const Foam::dictionary* Foam::dictionaryEntry::dictPtr() const
{
    return this;
}


Foam::dictionary* Foam::dictionaryEntry::dictPtr()
{
    return this;
}


const Foam::dictionary& Foam::dictionaryEntry::dict() const
{
    return *this;
}


Foam::dictionary& Foam::dictionaryEntry::dict()
{
    return *this;
}


// ************************************************************************* //
