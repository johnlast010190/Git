/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : dev
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
    (c) 1991-2005 OpenCFD Ltd. ICE Stroemungsfoschungs GmbH

Application
    funkySetFields

Description

Contributors/Copyright:
    2010, 2012-2014, 2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "PatchValueExpressionDriver.H"
#include "db/Time/timeSelector.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "helpers/printSwakVersion.H"
#include "repositories/RepositoryBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{

    Foam::timeSelector::addOptions(false);

#include "include/addRegionOption.H"

#include "appSnipplets/addLoadFunctionPlugins.H"

    argList::validOptions.insert("dict","<dictionary to use>");
    argList::validOptions.insert("cacheFields","");

#include "include/setRootCase.H"

    printSwakVersion();

    word dictName="funkySetBoundaryDict";
    if (args.options().found("dict")) {
        dictName=args.options()["dict"];
    }

    bool cacheFields=args.options().found("cacheFields");
    if (cacheFields) {
        WarningIn("main()")
            << "The current behaviour is to cache fields that were read from disc. "
                << "This may lead to unexpected behaviour as previous modifications "
                << "will not be visible."
                << endl;
            }

#include "include/createTime.H"
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);

#include "include/createNamedMesh.H"

#include "appSnipplets/loadFunctionPlugins.H"

    IOdictionary funkyDict
        (
            IOobject
            (
                dictName,
                runTime.system(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Foam::Info<< "Time = " << runTime.timeName() << Foam::endl;

        RepositoryBase::updateRepos();

        mesh.readUpdate();

        forAllIter(dictionary,funkyDict,it) {
            const dictionary &part=(*it).dict();

            word fieldName=part["field"];

            Info<< "\n\nPart: " << (*it).keyword()
                << " working on field " << fieldName << endl;

            IOdictionary field(
                IOobject
                (
                    fieldName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                )
            );
            // deregister the dictionary so that the field can work on itself
            field.checkOut();
            {
                // this way it doesn't matter that the file is not of the right class
                IFstream inStream(field.filePath());
                field.readHeader(inStream);
                field.readData(inStream);
            }

            List<dictionary> expressions(part.lookup("expressions"));

            forAll(expressions,expressionI) {
                const dictionary &expression=expressions[expressionI];

                word target(expression["target"]);
                word patchName(expression["patchName"]);
                exprString expr(
                    expression["expression"],
                    expression
                );
                Info<< "Setting " << target << " on " << patchName
                    << " the expression " << expr << endl;

                PatchValueExpressionDriver driver(expression,mesh);
                driver.setSearchBehaviour(
                    cacheFields,
                    false,
                    true             // search on disc
                );

                driver.clearVariables();
                driver.parse(expr);

                dictionary &patchDict=field.subDict("boundaryField").subDict(patchName);

                if (patchDict.found(target)) {
                    // Does not work (memory problem)
                    //                    patchDict.changeKeyword(keyType(target),keyType(word(target+"Old")),true);
                    if (patchDict.found(target+"Old")) {
                        patchDict.remove(target+"Old");
                    }
                    patchDict.changeKeyword(keyType(target),keyType(word(target+"Old")));
                }
                OStringStream result;
                string newEntry=driver.outputEntry();
                patchDict.set(target,newEntry.c_str());
            }

            {
                // this way the class is not overwritten
                word actualClass=field.headerClassName();

                OStringStream headerStream;
                field.writeHeader(headerStream);
                string newHeader=headerStream.str().replace("dictionary",actualClass);

                IFstream inStream(field.filePath());
                OFstream outStream(
                    field.filePath(),
#ifdef FOAM_DEV
            std::ios_base::out,
#endif
                    inStream.format(),
                    inStream.version(),
                    inStream.compression()
                );
                outStream << newHeader.c_str();
                field.writeData(outStream);
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
