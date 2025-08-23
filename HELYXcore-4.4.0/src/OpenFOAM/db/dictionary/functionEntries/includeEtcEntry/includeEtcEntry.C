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
    (c) 2017 Engys Ltd.
    (c) 2015-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "db/dictionary/functionEntries/includeEtcEntry/includeEtcEntry.H"
#include "global/etcFiles/etcFiles.H"
#include "primitives/strings/stringOps/stringOps.H"
#include "db/IOobject/IOobject.H"
#include "global/fileOperations/fileOperation/fileOperation.H"
#include "db/runTimeSelection/memberFunctions/addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::functionEntries::includeEtcEntry::log(false);


namespace Foam
{
namespace functionEntries
{
    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        includeEtcEntry,
        execute,
        dictionaryIstream,
        includeEtc
    );

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        includeEtcEntry,
        execute,
        primitiveEntryIstream,
        includeEtc
    );
}
}

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

Foam::fileName Foam::functionEntries::includeEtcEntry::resolveFile
(
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

    // Search the etc directories for the file
    return Foam::findEtcFile(fName);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::includeEtcEntry::mergeFile
(
    dictionary& parentDict,
    fileName rawFName
)
{

    const fileName fName(resolveFile(rawFName, parentDict));

    //IFstream ifs(fName);
    autoPtr<ISstream> ifsPtr(fileHandler().NewIFstream(fName));
    ISstream& ifs = ifsPtr();

    if (ifs)
    {
        if (Foam::functionEntries::includeEtcEntry::log)
        {
            Info<< fName << endl;
        }
        parentDict.read(ifs);
        return true;
    }
    else
    {
        FatalIOErrorInFunction
        (
            ifs
        )   << "Cannot open include file "
            << (ifs.name().size() ? ifs.name() : rawFName)
            << " while reading dictionary " << parentDict.name()
            << exit(FatalIOError);

        return false;
    }
}


bool Foam::functionEntries::includeEtcEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    const fileName rawName(is);
    const fileName fName(resolveFile(rawName, parentDict));

    //IFstream ifs(fName);
    autoPtr<ISstream> ifsPtr(fileHandler().NewIFstream(fName));
    ISstream& ifs = ifsPtr();

    if (ifs)
    {
        if (Foam::functionEntries::includeEtcEntry::log)
        {
            Info<< fName << endl;
        }
        // Cache the FoamFile entry if present
        dictionary foamFileDict;
        if (parentDict.found("FoamFile"))
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
        )   << "Cannot open etc file "
            << (ifs.name().size() ? ifs.name() : rawName)
            << " while reading dictionary " << parentDict.name()
            << exit(FatalIOError);

        return false;
    }
}


bool Foam::functionEntries::includeEtcEntry::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    const fileName rawName(is);
    const fileName fName(resolveFile(rawName, parentDict));

    //IFstream ifs(fName);
    autoPtr<ISstream> ifsPtr(fileHandler().NewIFstream(fName));
    ISstream& ifs = ifsPtr();

    if (ifs)
    {
        if (Foam::functionEntries::includeEtcEntry::log)
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
        )   << "Cannot open etc file "
            << (ifs.name().size() ? ifs.name() : rawName)
            << " while reading dictionary " << parentDict.name()
            << exit(FatalIOError);

        return false;
    }
}


// ************************************************************************* //
