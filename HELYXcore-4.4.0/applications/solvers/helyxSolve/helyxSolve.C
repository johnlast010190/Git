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
    (c) 2019-2024 Engys Ltd.

Application
    helyxSolve.C

Description
    General single- or multi-region runner for solver objects

SeeAlso
    Foam::helyxSolve

\*---------------------------------------------------------------------------*/

#include "helyxSolve.H"
#include "global/argList/argList.H"
#include "db/Time/timeSelector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
    // For Windows only
    Foam::argList::resetMaxOpenFileStreams();

    // Add postProcessing option
    Foam::argList::addBoolOption
    (
        argList::postProcessOptionName,
        "Execute functionObjects only"
    );
    Foam::argList::addOption
    (
        "instances",
        "list",
        "specify a list of instances to run"
    );
    if (argList::postProcess(argc, argv))
    {
        Foam::timeSelector::addOptions();
        #include "include/addRegionOption.H"
        #include "include/addFunctionObjectOptions.H"

        // Set functionObject post-processing mode
        functionObject::postProcess = true;
    }

    // Add profiling option
    #if !defined( WIN32 ) && !defined( WIN64 )
    #include "include/addProfilingOption.H"
    #endif

    // Set root case
    #include "include/setRootCase.H"

    fileName exeName(args.executable());

    Info<< "Creating time\n" << endl;
    Time time(Time::controlDictName, args);

    return helyxSolve(time).runSolver(args, exeName.nameLessExt());
}

// ************************************************************************* //
