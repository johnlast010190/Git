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
    solidDisplacementFoam

Group
    grpStressAnalysisSolvers

Description
    Transient segregated finite-volume solver of linear-elastic,
    small-strain deformation of a solid body, with optional thermal
    diffusion and thermal stresses.

    Simple linear elasticity structural analysis code.
    Solves for the displacement vector field D, also generating the
    stress tensor field sigma.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "solidDisplacementThermo/solidDisplacementThermo.H"
#include "primitives/bools/Switch/Switch.H"
#include  "cfdTools/general/fvOptions/fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define CREATE_CONTROL createControl.H
    #include "db/functionObjects/functionObjectList/postProcess.H"
    #undef CREATE_CONTROL

    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "include/createMesh.H"
    #include "createControls.H"
    #include "createFields.H"

    fv::options& fvOptions(fv::options::New(mesh));

    tmp<volScalarField> trho = thermo.rho();
    const volScalarField& rho = trho();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating displacement field\n" << endl;

    while (runTime.loop())
    {
        Info<< "Iteration: " << runTime.value() << nl << endl;

        #include "readSolidDisplacementFoamControls.H"

        int iCorr = 0;
        scalar initialResidual = 0;

        do
        {
            if (thermo.thermalStress())
            {
                volScalarField& T = Tptr();
                solve
                (
                    fvm::ddt(rho, Cp, T) == fvm::laplacian(kappa, T)
                  + fvOptions(rho*Cp, T)
                );
            }

            {
                fvVectorMatrix DEqn
                (
                    fvm::d2dt2(rho, D)
                 ==
                    fvm::laplacian(2*mu + lambda, D, "laplacian(DD,D)")
                  + divSigmaExp
                  + rho*fvOptions.d2dt2(D)
                );

                if (thermo.thermalStress())
                {
                    const volScalarField& T = Tptr();
                    DEqn += fvc::grad(threeKalpha*T);
                }

                initialResidual = DEqn.solve().max().initialResidual();

                if (!compactNormalStress)
                {
                    divSigmaExp = fvc::div(DEqn.flux());
                }
            }

            {
                volTensorField gradD(fvc::grad(D));
                sigmaD = mu*twoSymm(gradD) + (lambda*I)*tr(gradD);

                if (compactNormalStress)
                {
                    divSigmaExp = fvc::div
                    (
                        sigmaD - (2*mu + lambda)*gradD,
                        "div(sigmaD)"
                    );
                }
                else
                {
                    divSigmaExp += fvc::div(sigmaD);
                }
            }

        } while (initialResidual > convergenceTolerance && ++iCorr < nCorr);

        #include "calculateStress.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
