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
    sprayEngineFoam

Group
    grpLagrangianSolvers grpMovingMeshSolvers

Description
    Transient solver for compressible, turbulent engine flow with a spray
    particle cloud.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "engineTime/engineTime.H"
#include "engineMesh/engineMesh/engineMesh.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "clouds/derived/basicSprayCloud/basicSprayCloud.H"
#include "rhoMulticomponentThermo/rhoMulticomponentThermo.H"
#include "combustionModel/combustionModel.H"
#include "solverObjects/radiation/radiationSolver.H"
#include "SLGThermo/SLGThermo.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "cfdTools/general/fvOptions/fvOptions.H"

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
    #include "readEngineTimeControls.H"
    #include "createFields.H"
    const volScalarField& T = thermo.T();
    const volScalarField& psi = thermo.psi();
    #include "cfdTools/compressible/createRhoUf.H"
    #include "cfdTools/compressible/compressibleCourantNo.H"
    #include "cfdTools/general/include/setInitialDeltaT.H"
    #include "cfdTools/general/include/initContinuityErrs.H"
    #include "startSummary.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readEngineTimeControls.H"
        #include "cfdTools/compressible/compressibleCourantNo.H"
        #include "cfdTools/general/include/setDeltaT.H"

        runTime++;

        Info<< "Crank angle = " << runTime.theta() << " CA-deg" << endl;

        mesh.moveMesh();

        parcels.evolve();

        #include "rhoEqn.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            #include "YEqn.H"
            #include "EEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }

            fvOptions.correct();
        }

        #include "logSummary.H"

        rho = thermo.rho();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
