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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "db/dictionary/functionEntries/includeEntry/includeEntry.H"
#include "primitives/strings/stringOps/stringOps.H"
#include "db/IOobject/IOobject.H"
#include "global/fileOperations/fileOperation/fileOperation.H"
#include "db/runTimeSelection/memberFunctions/addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::functionEntries::includeEntry::log(false);


namespace Foam
{
namespace functionEntries
{
    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        includeEntry,
        execute,
        dictionaryIstream,
        include
    );

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        includeEntry,
        execute,
        primitiveEntryIstream,
        include
    );
}
}

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

Foam::fileName Foam::functionEntries::includeEntry::resolveFile
(
    Istream& is,
    const dictionary& dict
)
{
    fileName fName(is);
    // Substitute dictionary and environment variables.
    // Allow empty substitutions.
    stringOps::inplaceExpand(fName, dict, true, true);

    if (fName.empty() || fName.isAbsolute())
    {
        return fName;
    }

    // Relative name
    return fileName(is.name()).path()/fName;
}


Foam::fileName Foam::functionEntries::includeEntry::resolveFile
(
    const fileName& dir,
    const fileName& f,
    const dictionary& dict
)
{
    fileName fName(f);
    // Substitute dictionary and environment variables.
    // Allow empty substitutions.
    stringOps::inplaceExpand(fName, dict, true, true);

    if (fName.empty() || fName.isAbsolute())
    {
        return fName;
    }

    // Relative name
    return dir/fName;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::includeEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    const fileName rawName(is);
    const fileName fName = resolveFile(is.name().path(), rawName, parentDict);
    autoPtr<ISstream> ifsPtr(fileHandler().NewIFstream(fName));
    ISstream& ifs = ifsPtr();

    if (ifs)
    {
        if (Foam::functionEntries::includeEntry::log)
        {
            Info<< fName << endl;
        }

        // Cache the FoamFile entry if present
        dictionary foamFileDict;
        if (parentDict.found(IOobject::foamFile))
        {
            foamFileDict = parentDict.subDict(IOobject::foamFile);
        }

        // Read and clear the FoamFile entry
        parentDict.read(ifs);

        // Reinstate original FoamFile entry
        if (foamFileDict.size() != 0)
        {
            dictionary parentDictTmp(parentDict);
            parentDict.clear();
            parentDict.add(IOobject::foamFile, foamFileDict);
            parentDict += parentDictTmp;
        }
        return true;
    }
    else
    {
        FatalIOErrorInFunction
        (
            is
        )   << "Cannot open include file "
            << (ifs.name().size() ? ifs.name() : rawName)
            << " while reading dictionary " << parentDict.name()
            << exit(FatalIOError);

        return false;
    }
}


bool Foam::functionEntries::includeEntry::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    const fileName rawName(is);
    const fileName fName = resolveFile(is.name().path(), rawName, parentDict);
    autoPtr<ISstream> ifsPtr(fileHandler().NewIFstream(fName));
    ISstream& ifs = ifsPtr();

    if (ifs)
    {
        if (Foam::functionEntries::includeEntry::log)
        {
            Info<< fName << endl;
        }

        entry.read(parentDict, ifs);

        return true;
    }
    else
    {
        FatalIOErrorInFunction
        (
            is
        )   << "Cannot open include file "
            << (ifs.name().size() ? ifs.name() : rawName)
            << " while reading dictionary " << parentDict.name()
            << exit(FatalIOError);

        return false;
    }
}

// ************************************************************************* //
