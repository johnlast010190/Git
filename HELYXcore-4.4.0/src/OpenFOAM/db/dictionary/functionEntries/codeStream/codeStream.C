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
    (c) 2021-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "db/dictionary/functionEntries/codeStream/codeStream.H"
#include "db/dynamicLibrary/dynamicCode/dynamicCode.H"
#include "db/dynamicLibrary/dynamicCode/dynamicCodeContext.H"
#include "db/dynamicLibrary/dlLibraryTable/dlLibraryTable.H"
#include "db/regIOobject/regIOobject.H"
#include "include/OSspecific.H"
#include "db/runTimeSelection/memberFunctions/addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineTypeNameAndDebug(codeStream, 1);

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        codeStream,
        execute,
        dictionaryIstream,
        codeStream
    );

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        codeStream,
        execute,
        primitiveEntryIstream,
        codeStream
    );
}
}

const Foam::word Foam::functionEntries::codeStream::codeTemplateC =
    "codeStreamTemplate.C";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionEntries::codeStream::masterOnlyRead
(
    const dictionary& dict
)
{
    const dictionary& topDict = dict.topDict();

    if (debug)
    {
        Pout<< "codeStream: dictionary: " << dict.name()
            << " master-only-reading:" << topDict.global()
            << endl;
    }

    return topDict.global();
}


Foam::functionEntries::codeStream::streamingFunctionType
Foam::functionEntries::codeStream::getFunction
(
    const dictionary& parentDict,
    const dictionary& codeDict
)
{
    // get code, codeInclude, codeOptions
    dynamicCodeContext context(codeDict);

    // codeName: codeStream + _<sha1>
    // codeDir : _<sha1>
    const std::string sha1Str(context.sha1().str(true));
    dynamicCode dynCode("codeStream" + sha1Str, sha1Str);

    // Load library if not already loaded
    // Version information is encoded in the libPath (encoded with the SHA1)
    #if defined( WIN32 ) || defined( WIN64 )
        fileName const& libPath = dynCode.libPath();
    #else
        const fileName libPath = dynCode.libPath();
    #endif

    // See if library is loaded
    void* lib = libs.findLibrary(libPath);

    if (!lib)
    {
        Info<< "Using #codeStream with " << libPath << endl;
    }

    // Nothing loaded.
    // Avoid compilation if possible by loading an existing library.
    if (!lib)
    {
        // Cached access to dl libs.
        // Guarantees clean up upon destruction of Time.
        if (libs.open(libPath, false))
        {
            lib = libs.findLibrary(libPath);
        }
        else
        {
            // Uncached opening of libPath. Do not complain if cannot be loaded.
            lib = dlOpen(libPath, false);
        }
    }

    // Create library if required
    if (!lib)
    {
        const bool create =
            Pstream::master()
         || (regIOobject::fileModificationSkew <= 0);   // not NFS

        if (create)
        {
            if (!dynCode.upToDate(context))
            {
                // Filter with this context
                dynCode.reset(context);

                // Compile filtered C template
                dynCode.addCompileFile(codeTemplateC);

                // Define Make/options
                dynCode.setOptions(context.options());
                dynCode.setLinkLibs(context.libs());

                if (!dynCode.copyOrCreateFiles(true))
                {
                    FatalIOErrorInFunction
                    (
                        parentDict
                    )   << "Failed writing files for" << nl
                        << dynCode.libRelPath() << nl
                        << exit(FatalIOError);
                }
            }

            if (!dynCode.makeLibso())
            {
                FatalIOErrorInFunction
                (
                    parentDict
                )   << "Failed make " << dynCode.libRelPath() << nl
                    << exit(FatalIOError);
            }
        }

        // Only block if not master only reading of a global dictionary
        if
        (
           !masterOnlyRead(parentDict)
         && regIOobject::fileModificationSkew > 0
        )
        {
            // Since the library has only been compiled on the master the
            // other nodes need to pick this library up through NFS
            // We do this by just polling a few times using the
            // fileModificationSkew.
        #if defined( WIN32 ) || defined( WIN64 )
            off64_t mySize = fileSize(libPath,false);
        #else
            off64_t mySize = fileSize(libPath);
        #endif
            off64_t masterSize = mySize;
            Pstream::scatter(masterSize);

            if (debug)
            {
                Pout<< endl<< "on processor " << Pstream::myProcNo()
                    << " have masterSize:" << masterSize
                    << " and localSize:" << mySize
                    << endl;
            }


            if (mySize < masterSize)
            {
                if (debug)
                {
                    Pout<< "Local file " << libPath
                        << " not of same size (" << mySize
                        << ") as master ("
                        << masterSize << "). Waiting for "
                        << regIOobject::fileModificationSkew
                        << " seconds." << endl;
                }
                Foam::sleep(regIOobject::fileModificationSkew);

                // Recheck local size

            #if defined(WIN32) || defined(WIN64)
                mySize = Foam::fileSize(libPath,false);
            #else
                mySize = Foam::fileSize(libPath);
            #endif

                if (mySize < masterSize)
                {
                    FatalIOErrorInFunction
                    (
                        parentDict
                    )   << "Cannot read (NFS mounted) library " << nl
                        << libPath << nl
                        << "on processor " << Pstream::myProcNo()
                        << " detected size " << mySize
                        << " whereas master size is " << masterSize
                        << " bytes." << nl
                        << "If your case is not NFS mounted"
                        << " (so distributed) set fileModificationSkew"
                        << " to 0"
                        << exit(FatalIOError);
                }
            }

            if (debug)
            {
                Pout<< endl<< "on processor " << Pstream::myProcNo()
                    << " after waiting: have masterSize:" << masterSize
                    << " and localSize:" << mySize
                    << endl;
            }
        }

        if (libs.open(libPath, false))
        {
            if (debug)
            {
                Pout<< "Opening cached dictionary:" << libPath << endl;
            }

            lib = libs.findLibrary(libPath);
        }
        else
        {
            // Uncached opening of libPath
            if (debug)
            {
                Pout<< "Opening uncached dictionary:" << libPath << endl;
            }

            lib = dlOpen(libPath, true);
        }
    }

    if (!lib)
    {
        FatalIOErrorInFunction(parentDict)
            << "Failed loading library " << libPath << nl
            << "Did you add all libraries to the 'libs' entry"
            << " in system/controlDict?"
            << exit(FatalIOError);
    }

    bool haveLib = lib;
    if (!masterOnlyRead(parentDict))
    {
        reduce(haveLib, andOp<bool>());
    }

    if (!haveLib)
    {
        FatalIOErrorInFunction
        (
            parentDict
        )   << "Failed loading library " << libPath
            << " on some processors."
            << exit(FatalIOError);
    }

    // Find the function handle in the library
    const streamingFunctionType function =
        reinterpret_cast<streamingFunctionType>
        (
            dlSym(lib, dynCode.codeName())
        );

    if (!function)
    {
        FatalIOErrorInFunction
        (
            parentDict
        )   << "Failed looking up symbol " << dynCode.codeName()
            << " in library " << lib << exit(FatalIOError);
    }

    return function;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::codeStream::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    Info<< "Using #codeStream at line " << is.lineNumber()
        << " in file " <<  parentDict.name() << endl;

    dynamicCode::checkSecurity
    (
        "functionEntries::codeStream::execute(..)",
        parentDict
    );

    // Get code dictionary.
    // Must reference parent for stringOps::expand to work nicely.
    const dictionary codeDict("#codeStream", parentDict, is);

    const streamingFunctionType function = getFunction(parentDict, codeDict);

    // Use function to write stream
    OStringStream os(is.format());
    (*function)(os, parentDict);

    // Get the entry from this stream
    IStringStream resultStream(os.str());
    entry.read(parentDict, resultStream);

    return true;
}


bool Foam::functionEntries::codeStream::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    Info<< "Using #codeStream at line " << is.lineNumber()
        << " in file " <<  parentDict.name() << endl;

    dynamicCode::checkSecurity
    (
        "functionEntries::codeStream::execute(..)",
        parentDict
    );

    // Get code dictionary.
    // Must reference parent for stringOps::expand to work nicely.
    const dictionary codeDict("#codeStream", parentDict, is);

    const streamingFunctionType function = getFunction(parentDict, codeDict);

    // Use function to write stream
    OStringStream os(is.format());
    (*function)(os, parentDict);

    // Get the entry from this stream
    IStringStream resultStream(os.str());
    parentDict.read(resultStream);

    return true;
}


// ************************************************************************* //
