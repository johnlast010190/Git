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
    (c) 2011-2016 OpenFOAM Foundation

Application
    simpleFoam

Description
    Executes solverObject-based fvOptions. Incompressible flow solvers
    with turbulence modelling.


\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "cfdTools/general/solutionControl/simpleControl/simpleControl.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "cfdTools/general/fvOptions/fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "include/setRootCase.H"

    DeprecationWarningInFunction
    (
        "solverObjectFoam",
        "solver",
        40300,
        "Please use the helyxSolve application instead."
    );

    #include "include/createTime.H"
    #include "include/createMesh.H"
    #include "cfdTools/general/include/createSimpleControl.H"
    #include "cfdTools/general/include/createPimpleControl.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    if (solutionMode == "SIMPLE" || solutionMode == "Simple" || solutionMode == "simple")
    {
        while (simple.loop())
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;

            fvOptions.correct();

            runTime.write();

            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;
        }
    }
    else if (solutionMode == "PIMPLE" || solutionMode == "Pimple" || solutionMode == "pimple")
    {
        #include "cfdTools/general/include/createTimeControls.H"
        #include "cfdTools/incompressible/CourantNo.H"
        #include "cfdTools/general/include/setInitialDeltaT.H"

        while (runTime.run())
        {
            #include "cfdTools/general/include/readTimeControls.H"
            #include "cfdTools/incompressible/CourantNo.H"
            #include "cfdTools/general/include/setDeltaT.H"

            runTime++;

            Info<< "Time = " << runTime.timeName() << nl << endl;

            fvOptions.correct();

            runTime.write();

            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
