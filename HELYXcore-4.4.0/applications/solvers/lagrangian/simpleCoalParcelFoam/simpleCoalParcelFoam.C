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
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2025 Engys Ltd.

Application
    simpleReactingParcelFoam

Group
    grpLagrangianSolvers

Description
    Steady state solver for compressible, turbulent flow with coal particle
    clouds and optional sources/constraints.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "coalCloud/coalCloud.H"
#include "rhoMulticomponentThermo/rhoMulticomponentThermo.H"
#include "combustionModel/combustionModel.H"
#include "solverObjects/radiation/radiationSolver.H"
#include "cfdTools/general/porosityModel/porosityModel/IOporosityModelList.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "SLGThermo/SLGThermo.H"
#include "cfdTools/general/solutionControl/simpleControl/simpleControl.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "db/functionObjects/functionObjectList/postProcess.H"

    #include "include/setRootCase.H"
    DeprecationWarningInFunction
    (
        "simpleCoalParcelFoam",
        "solver",
        40400,
        "Please use the helyxSolve application instead."
    );
    #include "include/createTime.H"
    #include "include/createMesh.H"
    #include "cfdTools/general/solutionControl/createControl.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "cfdTools/general/include/initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        parcels.evolve();

        // --- Pressure-velocity SIMPLE corrector loop
        {
            #include "UEqn.H"
            #include "YEqn.H"
            #include "EEqn.H"
            #include "pEqn.H"
        }

        turbulence->correct();

        fvOptions.correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
