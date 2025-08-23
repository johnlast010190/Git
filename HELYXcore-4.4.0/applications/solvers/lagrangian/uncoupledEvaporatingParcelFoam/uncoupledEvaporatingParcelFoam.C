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
    (c) 2010-2016 Engys Ltd.
    (c) 2011-2020 OpenFOAM Foundation

Application
    reactingParcelFoam

Group
    grpLagrangianSolvers

Description
    Transient PIMPLE solver for compressible, laminar or turbulent flow with
    reacting multiphase Lagrangian parcels, including run-time selectable
    finite volume options, e.g. sources, constraints

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "clouds/derived/basicReactingMultiphaseCloud/basicReactingMultiphaseCloud.H"
#include "rhoMulticomponentThermo/rhoMulticomponentThermo.H"
#include "combustionModel/combustionModel.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "SLGThermo/SLGThermo.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "include/setRootCase.H"

    DeprecationWarningInFunction
    (
        "uncoupledEvaporatingParcelFoam",
        "solver",
        40400,
        "Please use the helyxSolve application instead."
    );

    #include "include/createTime.H"
    #include "include/createMesh.H"
    #include "cfdTools/general/include/readGravitationalAcceleration.H"

    pimpleControl pimple(mesh);

    #include "createFields.H"
    #include "createClouds.H"
    #include "cfdTools/general/include/createFvOptions.H"
    #include "cfdTools/general/include/createTimeControls.H"
    #include "cfdTools/compressible/compressibleCourantNo.H"
    #include "cfdTools/general/include/setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "cfdTools/general/include/readTimeControls.H"
        #include "cfdTools/compressible/compressibleCourantNo.H"
        #include "cfdTools/general/include/setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        parcels.evolve();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
