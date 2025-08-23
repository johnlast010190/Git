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
    reactingParcelFoam

Group
    grpLagrangianSolvers

Description
    Transient solver for compressible, turbulent flow with a reacting,
    multiphase particle cloud, and surface film modelling.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "clouds/derived/basicReactingMultiphaseCloud/basicReactingMultiphaseCloud.H"
#include "surfaceFilmModels/surfaceFilmModel/surfaceFilmModel.H"
#include "rhoMulticomponentThermo/rhoMulticomponentThermo.H"
#include "combustionModel/combustionModel.H"
#include "solverObjects/radiation/radiationSolver.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "SLGThermo/SLGThermo.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "cfdTools/general/pressureControl/pressureControl.H"
#include "finiteVolume/ddtSchemes/localEulerDdtScheme/localEulerDdtScheme.H"
#include "finiteVolume/fvc/fvcSmooth/fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "db/functionObjects/functionObjectList/postProcess.H"

    #include "include/setRootCase.H"

    DeprecationWarningInFunction
    (
        "reactingParcelFoam",
        "solver",
        40300,
        "Please use the helyxSolve application instead."
    );

    #include "include/createTime.H"
    #include "include/createMesh.H"
    #include "cfdTools/general/solutionControl/createControl.H"
    #include "cfdTools/general/include/createTimeControls.H"
    #include "createFields.H"
    Switch solvePrimaryRegion
    (
        pimple.dict().lookupOrDefault<Switch>("solvePrimaryRegion", true)
    );
    const volScalarField& T = thermo.T();
    const volScalarField& psi = thermo.psi();
    regionModels::surfaceFilmModel& surfaceFilm = tsurfaceFilm();
    #include "cfdTools/general/include/initContinuityErrs.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "cfdTools/compressible/compressibleCourantNo.H"
        #include "cfdTools/general/include/setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "cfdTools/general/include/readTimeControls.H"
        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "cfdTools/compressible/compressibleCourantNo.H"
            if (adjustTimeStep)
            {
                const scalar maxDeltaTFact =
                    min
                    (
                        maxCo/(CoNum + SMALL),
                        maxCo/(surfaceFilm.CourantNumber() + SMALL)
                    );
                const scalar deltaTFact =
                    min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);
                runTime.setDeltaT
                (
                    min(deltaTFact*runTime.deltaTValue(), maxDeltaT)
                );
                Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
            }
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        parcels.evolve();
        surfaceFilm.evolve();

        if (solvePrimaryRegion && !pimple.SIMPLErho())
        {
            #include "rhoEqn.H"
        }

        // --- PIMPLE loop
        while (solvePrimaryRegion && pimple.loop())
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
                thermo.correct();
            }
        }

        rho = thermo.rho();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
