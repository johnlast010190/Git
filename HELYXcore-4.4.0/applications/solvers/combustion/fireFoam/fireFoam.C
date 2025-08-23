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
    (c) 2025 Engys Ltd.

Application
    fireFoam

Group
    grpCombustionSolvers

Description
    Transient solver for fires and turbulent diffusion flames with reacting
    particle clouds, surface film and pyrolysis modelling.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "clouds/derived/basicReactingCloud/basicReactingCloud.H"
#include "surfaceFilmModels/surfaceFilmModel/surfaceFilmModel.H"
#include "pyrolysisModels/pyrolysisModel/pyrolysisModelCollection.H"
#include "solverObjects/radiation/radiationSolver.H"
#include "SLGThermo/SLGThermo.H"
#include "solidChemistryModel/solidChemistryModel.H"
#include "fluidMulticomponentThermo/fluidMulticomponentThermo.H"
#include "combustionModel/combustionModel.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "cfdTools/general/fvOptions/fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "db/functionObjects/functionObjectList/postProcess.H"

    #include "include/setRootCase.H"

    DeprecationWarningInFunction
    (
        "fireFoam",
        "solver",
        40300,
        "Please use the helyxSolve application instead."
    );

    #include "include/createTime.H"
    #include "include/createMesh.H"
    #include "cfdTools/general/solutionControl/createControl.H"
    #include "createFields.H"
    const volScalarField& psi = thermo.psi();
    const volScalarField& T = thermo.T();
    regionModels::surfaceFilmModel& surfaceFilm = tsurfaceFilm();
    #include "cfdTools/general/include/initContinuityErrs.H"
    #include "cfdTools/general/include/createTimeControls.H"
    #include "cfdTools/compressible/compressibleCourantNo.H"
    #include "cfdTools/general/include/setInitialDeltaT.H"
    scalar maxDi = pyrolysis.maxDiff();

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "cfdTools/general/include/readTimeControls.H"
        #include "cfdTools/compressible/compressibleCourantNo.H"
        scalar DiNum =  pyrolysis.solidRegionDiffNo();

        // Reset the timestep to maintain a constant maximum Courant numbers.
        // Reduction of time-step is immediate, but increase is damped to avoid
        // unstable oscillations.
        if (adjustTimeStep)
        {
            if (CoNum == -GREAT)
            {
                CoNum = SMALL;
            }
            const scalar TFactorFluid = maxCo/(CoNum + SMALL);
            const scalar TFactorSolid = maxDi/(DiNum + SMALL);
            const scalar TFactorFilm =
                maxCo/(surfaceFilm.CourantNumber() + SMALL);
            const scalar dt0 = runTime.deltaTValue();
            runTime.setDeltaT
            (
                min
                (
                    dt0
                   *min
                    (
                        min(TFactorFluid, min(TFactorFilm, TFactorSolid)),
                        1.2
                    ),
                    maxDeltaT
                )
            );

            // Exit Helyx if deltaT is below user secified minDeltaT
            if (runTime.deltaTValue() < minDeltaT)
            {
                FatalErrorInFunction
                    << "minDeltaT = " << minDeltaT << nl
                    << "The computed deltaT is lower than the minDeltaT"
                    << exit(FatalError);
            }

        }
        #include "cfdTools/general/include/setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        parcels.evolve();

        surfaceFilm.evolve();

        if (solvePyrolysisRegion)
        {
            pyrolysis.evolve();
        }

        if (solvePrimaryRegion)
        {
            #include "rhoEqn.H"

            // --- PIMPLE loop
            while (pimple.loop())
            {
                #include "UEqn.H"
                #include "YEEqn.H"

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

            rho = thermo.rho();
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
