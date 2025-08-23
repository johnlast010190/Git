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
    (c) 2013-2016 OpenFOAM Foundation

Application
    compressibleMultiphaseInterFoam

Group
    grpMultiphaseSolvers

Description
    Solver for n compressible, non-isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Turbulence modelling is generic, i.e.  laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "multiphaseMixtureThermo.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "db/functionObjects/functionObjectList/postProcess.H"

    #include "include/setRootCase.H"

    DeprecationWarningInFunction
    (
        "compressibleMultiphaseInterFoam",
        "solver",
        40400,
        "Please use the helyxSolve application instead."
    );

    #include "include/createTime.H"
    #include "include/createMesh.H"
    #include "cfdTools/general/solutionControl/createControl.H"
    #include "cfdTools/general/include/createTimeControls.H"
    #include "createFields.H"
    #include "cfdTools/incompressible/CourantNo.H"
    #include "cfdTools/general/include/setInitialDeltaT.H"

    volScalarField& p = mixture.p();
    volScalarField& T = mixture.T();

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "cfdTools/general/include/readTimeControls.H"
        #include "cfdTools/incompressible/CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            mixture.solve();

            solve(fvm::ddt(rho) + fvc::div(mixture.rhoPhi()));

            #include "UEqn.H"
            #include "TEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
