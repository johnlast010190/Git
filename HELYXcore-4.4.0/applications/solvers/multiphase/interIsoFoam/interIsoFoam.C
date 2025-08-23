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
    (c) 2017 OpenCFD Ltd.
    (c) 2016 DHI
    (c) 2011-2017 OpenFOAM Foundation

Application
    interIsoFoam

Group
    grpMultiphaseSolvers

Description
    Solver derived from interFoam for 2 incompressible, isothermal immiscible
    fluids using the iso-advector phase-fraction based interface capturing
    approach.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

    For a two-fluid approach see twoPhaseEulerFoam.

    Reference:
        \verbatim
            Roenby, J., Bredmose, H. and Jasak, H. (2016).
            A computational method for sharp interface advection
            Royal Society Open Science, 3
            doi 10.1098/rsos.160405
        \endverbatim

    isoAdvector code supplied by Johan Roenby, DHI (2016)

\*---------------------------------------------------------------------------*/

#include "fvMatrices/solvers/isoAdvection/isoAdvection/isoAdvection.H"
#include "cfdTools/general/include/fvCFD.H"
#include "algorithms/subCycle/subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "cfdTools/general/CorrectPhi/CorrectPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "db/functionObjects/functionObjectList/postProcess.H"

    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "include/createMesh.H"
    #include "cfdTools/general/solutionControl/createControl.H"
    #include "cfdTools/general/include/createTimeControls.H"
    #include "cfdTools/general/include/initContinuityErrs.H"
    #include "createFields.H"
    #include "correctPhi.H"

    turbulence->validate();

    #include "cfdTools/general/include/readTimeControls.H"
    #include "cfdTools/incompressible/CourantNo.H"
    #include "cfdTools/general/include/setInitialDeltaT.H"

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
            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct();

            if (pimple.frozenFlow())
            {
                continue;
            }

            #include "UEqn.H"

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

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
