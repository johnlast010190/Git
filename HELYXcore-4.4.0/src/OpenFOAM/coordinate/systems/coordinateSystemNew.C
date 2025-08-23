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
    (c) 2024 Engys Ltd.
    (c) 2018 OpenCFD Ltd.
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "db/objectRegistry/objectRegistry.H"
#include "coordinate/systems/cartesianCS.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::dictionary* Foam::coordinateSystem::subDictCompat
(
    const dictionary* dictPtr
)
{
    if (dictPtr)
    {
        // Non-recursive, no pattern matching in the search
        const auto finder =
            dictPtr->csearch(coordinateSystem::typeName_(), keyType::LITERAL);

        if (finder.isDict())
        {
            return finder.dictPtr();
        }
        else if (finder.found())
        {
            const word csName(finder.ref().stream());

            // Deprecated, unsupported syntax

            DeprecationIOWarningInFunction
            (
                *dictPtr, "coordinateSystem", "keyword", 30000, ""
            );
        }
    }

    return dictPtr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::coordinateSystem> Foam::coordinateSystem::New
(
    word modelType,
    const objectRegistry& obr,
    const dictionary& dict
)
{
    if (modelType.empty())
    {
        modelType = coordSystem::cartesian::typeName_();
    }

    auto cstrIter1 = registryConstructorTable_().find(modelType);

    if (cstrIter1 != registryConstructorTable_().end())
    {
        return autoPtr<coordinateSystem>(cstrIter1->second(obr, dict));
    }

    const auto ctor =
        ctorTableLookup
        (
            "coordinate system type",
            dictionaryConstructorTable_(),
            modelType
        );
    return autoPtr<coordinateSystem>(ctor(dict));
}


Foam::autoPtr<Foam::coordinateSystem> Foam::coordinateSystem::New
(
    word modelType,
    const dictionary& dict
)
{
    if (modelType.empty())
    {
        modelType = coordSystem::cartesian::typeName_();
    }

    const auto ctor =
        ctorTableLookup
        (
            "coordinate system type",
            dictionaryConstructorTable_(),
            modelType
        );
    return autoPtr<coordinateSystem>(ctor(dict));
}


Foam::autoPtr<Foam::coordinateSystem> Foam::coordinateSystem::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& dictName
)
{
    const dictionary* dictPtr = &dict;

    if (dictName.size())
    {
        dictPtr = &(dictPtr->subDict(dictName));
    }
    else
    {
        // Use 'coordinateSystem' subDict if present
        dictPtr = coordinateSystem::subDictCompat(dictPtr);
    }

    word modelType = dictPtr->lookupOrDefault<word>
    (
        "type",
        coordSystem::cartesian::typeName_()
    );

    return coordinateSystem::New(modelType, obr, *dictPtr);
}


Foam::autoPtr<Foam::coordinateSystem> Foam::coordinateSystem::New
(
    const dictionary& dict,
    const word& dictName
)
{
    const dictionary* dictPtr = &dict;

    if (dictName.size())
    {
        dictPtr = &(dictPtr->subDict(dictName));
    }
    else
    {
        // Use 'coordinateSystem' subDict if present
        dictPtr = coordinateSystem::subDictCompat(dictPtr);
    }

    word modelType = dictPtr->lookupOrDefault<word>
    (
        "type",
        coordSystem::cartesian::typeName_()
    );

    return coordinateSystem::New(modelType, *dictPtr);
}


Foam::autoPtr<Foam::coordinateSystem> Foam::coordinateSystem::New(Istream& is)
{
    const word csName(is);
    const dictionary dict(is);

    auto cs = coordinateSystem::New(dict, word::null);
    cs->rename(csName);

    return cs;
}


// ************************************************************************* //
