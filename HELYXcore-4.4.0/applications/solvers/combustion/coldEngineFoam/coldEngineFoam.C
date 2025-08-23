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

Application
    coldEngineFoam

Group
    grpCombustionSolvers grpMovingMeshSolvers

Description
    Solver for cold-flow in internal combustion engines.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "engineTime/engineTime.H"
#include "engineMesh/engineMesh/engineMesh.H"
#include "fluidThermo/fluidThermo.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define CREATE_TIME include/createEngineTime.H
    #define CREATE_MESH include/createEngineMesh.H
    #include "db/functionObjects/functionObjectList/postProcess.H"

    #include "include/setRootCase.H"
    #include "include/createEngineTime.H"
    #include "include/createEngineMesh.H"
    #include "cfdTools/general/solutionControl/createControl.H"
    #include "createFields.H"

    const volScalarField& psi = thermo.psi();
    const volScalarField& T = thermo.T();

    #include "cfdTools/compressible/createRhoUf.H"
    #include "cfdTools/general/include/initContinuityErrs.H"
    #include "engineFoam/readEngineTimeControls.H"
    #include "cfdTools/compressible/compressibleCourantNo.H"
    #include "cfdTools/general/include/setInitialDeltaT.H"
    #include "startSummary.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "engineFoam/readEngineTimeControls.H"
        #include "cfdTools/compressible/compressibleCourantNo.H"
        #include "cfdTools/general/include/setDeltaT.H"
        #include "readControls.H"

        runTime++;

        Info<< "Crank angle = " << runTime.theta() << " CA-deg"
            << endl;

        mesh.moveMesh();

        #include "cfdTools/compressible/rhoEqn.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "engineFoam/UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "rhoPimpleFoam/rhoPimpleFoam/EEqn.H"
                #include "engineFoam/pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }

            fvOptions.correct();
        }

        runTime.write();

        #include "logSummary.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
