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
    (c) 2011-2023 OpenFOAM Foundation

Application
    chemFoam

Group
    grpCombustionSolvers

Description
    Solver for chemistry problems, designed for use on single cell cases to
    provide comparison against other chemistry solvers, that uses a single cell
    mesh, and fields created from the initial conditions.


\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "fluidMulticomponentThermo/fluidMulticomponentThermo.H"
#include "basicChemistryModel/basicChemistryModel.H"
#include "mixtures/multicomponentMixture/multicomponentMixture.H"
#include "chemistrySolver/chemistrySolver/chemistrySolver.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "mixtures/speciesMassFractions/speciesMassFractions.H"
#include "meshes/meshShapes/cellModeller/cellModeller.H"
#include "materialModels/baseModels/materialModels.H"
#include "materialModels/materialTables/materialTables.H"

scalarList W(const fluidMulticomponentThermo& thermo)
{
    const basicSpecieMixture& composition = thermo.composition();
    scalarList W(composition.Y().size());
    forAll(W, i)
    {
        W[i] = composition.W(i);
    }
    return W;
}


scalar h0
(
    const fluidMulticomponentThermo& thermo,
    const scalarList& Y,
    const scalar p,
    const scalar T
)
{
    const basicSpecieMixture& composition = thermo.composition();
    scalar h0 = 0;
    if (basicThermo::dictName == basicThermo::matDictName)
    {
        forAll(Y, i)
        {
            h0 +=
                Y[i]
               *thermo.materials()
                (
                    hsModel::typeName,
                    word::null,
                    thermo.composition().Y()[i].name()
                )[0];
        }
    }
    else
    {
        forAll(Y, i)
        {
            h0 += Y[i]*composition.hs(i, p, T);
        }
    }
    return h0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();

    #define CREATE_MESH createSingleCellMesh.H
    #define NO_CONTROL
    #include "db/functionObjects/functionObjectList/postProcess.H"

    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "createSingleCellMesh.H"
    #include "createFields.H"

    basicChemistryModel& chemistry = pChemistry();
    scalar dtChem = min(chemistry.deltaTChem()[0], runTime.deltaT().value());
    speciesMassFractions& composition = thermo.composition();
    PtrList<volScalarField>& Y = composition.Y();
    volScalarField& p = thermo.p();

    #include "readInitialConditions.H"

    Switch adjustTimeStep(runTime.controlDict().lookup("adjustTimeStep"));
    scalar maxDeltaT, minDeltaT;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        adjustTimeStep = runTime.controlDict().lookup<Switch>("adjustTimeStep");
        maxDeltaT = runTime.controlDict().lookup<scalar>("maxDeltaT");
        minDeltaT =
            runTime.controlDict().lookupOrDefault<scalar>("minDeltaT", SMALL);

        if (adjustTimeStep)
        {
            runTime.setDeltaT(min(dtChem, maxDeltaT));
            Info<< "deltaT = " <<  runTime.deltaT().value() << endl;

            // exit Helyx if deltaT is below user secified minDeltaT
            if (runTime.deltaTValue() < minDeltaT)
            {
                FatalErrorInFunction
                    << "minDeltaT = " << minDeltaT << nl
                    << "The computed deltaT is lower than the minDeltaT"
                    << exit(FatalError);
            }
        }

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Solve chemistry
        dtChem = chemistry.solve(runTime.deltaT().value());
        scalar Qdot = chemistry.Qdot()()[0]/rho[0];
        integratedHeat += Qdot*runTime.deltaT().value();

        // Species equations
        {
            forAll(Y, specieI)
            {
                volScalarField& Yi = Y[specieI];
                solve
                (
                    fvm::ddt(rho, Yi) - chemistry.RR()[specieI],
                    mesh.solution().solver("Yi")
                );
            }
        }

        // Enthalpy equation
        {
            volScalarField& h = thermo.he();
            if (constProp == "volume")
            {
                h[0] = u0 + p[0]/rho[0] + integratedHeat;
            }
            else
            {
                h[0] = h0 + integratedHeat;
            }
            thermo.correct();
        }

        // Pressure equation
        {
            rho = thermo.rho();
            if (constProp == "volume")
            {
                scalar invW = 0.0;
                forAll(Y, i)
                {
                    invW += Y[i][0]/W[i];
                }
                Rspecific[0] = 1000.0*constant::physicoChemical::R.value()*invW;
                p[0] = rho0*Rspecific[0]*(thermo.T()[0] + thermo.TRefValue());
                rho[0] = rho0;
            }
        }

        runTime.write();

        Info<< "Qdot = " << Qdot
            << ", T = " << thermo.T()[0]
            << ", p = " << thermo.p()[0]
            << ", " << Y[0].name() << " = " << Y[0][0]
            << endl;

        post<< runTime.value() << token::TAB << thermo.T()[0] << token::TAB
            << thermo.p()[0] << endl;

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "Number of steps = " << runTime.timeIndex() << nl
        << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
