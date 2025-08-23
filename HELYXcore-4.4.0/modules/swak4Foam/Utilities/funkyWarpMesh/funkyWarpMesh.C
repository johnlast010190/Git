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
    (c) 2024 Engys Ltd.

Application
    funkyWarpMesh

Description
    Warp mesh with new coordinates calculated for the points

Contributors/Copyright:
    2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "FieldValueExpressionDriver.H"
#include "helpers/printSwakVersion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
#include "include/addRegionOption.H"

#include "appSnipplets/addLoadFunctionPlugins.H"

    argList::validOptions.insert("overwrite", "");
    argList::validOptions.insert("expression","expression to write");
    argList::validOptions.insert("dictExt","Extension to the dictionary");
    argList::validOptions.insert("relative", "");
    argList::validOptions.insert("allowFunctionObjects","");

#include "include/setRootCase.H"

    printSwakVersion();

#include "include/createTime.H"

    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);

#include "include/createNamedMesh.H"

#include "appSnipplets/loadFunctionPlugins.H"

    if (!args.options().found("allowFunctionObjects")) {
        runTime.functionObjects().off();
    } else {
        runTime.functionObjects().start();
    }

    const word oldInstance = mesh.pointsInstance();

    bool overwrite    = args.optionFound("overwrite");
    bool relative     = false;

    pointField newPoints;

    if (!overwrite)
    {
        runTime++;
    }

    if (args.optionFound("expression")) {
        if (args.optionFound("dictExt")) {
            FatalErrorIn(args.executable())
                << "Can't specify 'dictExt' and 'expression' at the same time"
                    << endl
                    << exit(FatalError);

        }
        relative=args.optionFound("relative");
        exprString expression(
            args.options()["expression"],
            dictionary::null
        );
        FieldValueExpressionDriver driver(
            runTime.timeName(),
            runTime,
            mesh
        );
        // no clearVariables needed here
        driver.parse(expression);
        if (!driver.resultIsTyp<pointVectorField>()) {
            FatalErrorIn(args.executable())
                << "Expression " << expression
                    << " does not evaluate to a pointVectorField but a "
                    << driver.typ()
                    << endl
                    << exit(FatalError);
        }
        newPoints=driver.getResult<pointVectorField>().internalField();
    } else {
        Info<< "Dictionary mode" << nl << endl;
        if (args.optionFound("relative")) {
            FatalErrorIn(args.executable())
                << "Option 'relative' not allowed in dictionary-mode"
                    << endl
                    << exit(FatalError);
        }
        word dictName("funkyWarpMeshDict");
        if (args.optionFound("dictExt")) {
            dictName+="."+word(args.options()["dictExt"]);
        }
        IOdictionary warpDict
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
        const word mode(warpDict.lookup("mode"));
        if (mode=="set") {
            relative=warpDict.lookup<bool>("relative");
            exprString expression(
                warpDict.lookup("expression"),
                warpDict
            );
            FieldValueExpressionDriver driver(
                runTime.timeName(),
                runTime,
                mesh
            );
            driver.readVariablesAndTables(warpDict);
            driver.clearVariables();

            driver.parse(expression);
            if (!driver.resultIsTyp<pointVectorField>()) {
                FatalErrorIn(args.executable())
                    << "Expression " << expression
                        << " does not evaluate to a pointVectorField but a "
                        << driver.typ()
                        << endl
                        << exit(FatalError);
            }
            newPoints=driver.getResult<pointVectorField>().internalField();
        } else if (mode=="move") {
            notImplemented(args.executable()+" mode: move");
        } else {
            FatalErrorIn(args.executable())
                << "Possible values for 'mode' are 'set' or 'move'"
                    << endl
                    << exit(FatalError);
        }
    }

    if (relative) {
        newPoints += mesh.points();
    }

    mesh.movePoints(newPoints);

    // Write mesh
    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }
    Info<< nl << "Writing polyMesh to time " << runTime.timeName() << endl;

    IOstream::defaultPrecision(15);

    // Bypass runTime write (since only writes at outputTime)
    if
    (
       !runTime.objectRegistry::writeObject
        (
            runTime.writeFormat(),
            IOstream::currentVersion,
            runTime.writeCompression()
#ifdef FOAM_REGIOOBJECT_WRITEOBJECT_WITH_VALID
            ,true
#endif
        )
    )
    {
        FatalErrorIn(args.executable())
            << "Failed writing polyMesh."
            << exit(FatalError);
    }

    // Write points goes here

    // Write fields
    runTime.write();

    Info<< "End\n" << endl;

    Info<< nl << "Now run 'checkMesh' before you do anything else"
        << nl << endl;

    return 0;
}


// ************************************************************************* //
