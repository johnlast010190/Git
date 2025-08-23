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
    (c) 2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Application
    WatersKing

Description
    Analytical solution for the start-up planar Poiseuille flow of an
    Oldroyd-B fluid.

    Used in the planarPoiseuille example.

    References:
    \verbatim
        Waters, N. D., & King, M. J. (1970).
        Unsteady flow of an elasto-viscous liquid.
        Rheologica Acta, 9, 345-355.

        Amoreira, L. J., & Oliveira, P. J. (2010).
        Comparison of different formulations for the numerical
        calculation of unsteady incompressible viscoelastic fluid
        flow. Adv. Appl. Math. Mech, 4, 483-502.
    \endverbatim

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "turbulentTransportModels/turbulentTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "include/createMesh.H"
    #include "createFields.H"

    const scalar h = mesh.bounds().span().y();
    Info<< "Height from centreline to wall = " << h << endl;

    label centrelineID = mesh.boundary().findPatchID("centreline");
    const vector patchToCell =
        mesh.boundary()[centrelineID].Cf()[0]
      - mesh.C()[mesh.findNearestCell(location)];

    const scalar y = patchToCell.y()/h;
    Info<< "Normalised distance from centreline = " << y << nl << endl;

    const scalar nu0 = nu1 + nu2;
    const scalar E = lambda*nu0/(rho*sqr(h));
    const scalar beta = nu2/nu0;
    const scalar UInf = K*sqr(h)/3.0/nu0;

    Info<< "Waters and King parameters:" << nl
        << "E =    " << E << nl
        << "beta = " << beta << nl
        << "K =    " << K << nl
        << "UInf = " << UInf << nl << endl;

    label order = 8;

    scalarField ak(order, 0);
    scalarField bk(order, 0);
    scalarField ck(order, 0);
    scalarField B(order, 0);

    forAll(ak, i)
    {
        scalar k = i + 1;
        ak[i] = (2.0*k - 1)/2.0*constant::mathematical::pi*Foam::sqrt(E);
        bk[i] = (1.0 + beta*sqr(ak[i]))/2.0;
        ck[i] = Foam::sqrt(mag(sqr(bk[i]) - sqr(ak[i])));
        B[i]  = 48*Foam::pow(-1, k)
               /Foam::pow((2*k - 1)*constant::mathematical::pi, 3)
               *Foam::cos((2*k - 1)*constant::mathematical::pi*y/2);
    }

    scalarField A(order, 0);
    OFstream file(runTime.path()/"WatersKing.dat");
    const scalar LOGVGREAT = Foam::log(VGREAT);
    bool doResize = false;
    while (!runTime.end())
    {
        scalar t = runTime.timeOutputValue()/lambda;
        forAll(A, i)
        {
            if (bk[i]*t < LOGVGREAT)
            {
                if (bk[i] >= ak[i])
                {
                    A[i] = (bk[i] - sqr(ak[i]))/ck[i]*Foam::sinh(ck[i]*t)
                    + Foam::cosh(ck[i]*t);
                }
                else
                {
                    A[i] = (bk[i] - sqr(ak[i]))/ck[i]*Foam::sin(ck[i]*t)
                         + Foam::cos(ck[i]*t);
                }
                A[i] *= Foam::exp(-bk[i]*t);
            }
            else
            {
                Info<< "Coefficient A[" << order << "] = 0" << endl;
                order = i;
                doResize = true;
            }
        }

        if (doResize)
        {
            Info<< "Resizing A and B to " << order << endl;
            A.resize(order);
            B.resize(order);
            doResize = false;
        }

        scalar U = UInf*(1.5*(1 - sqr(y)) + sum(A*B));
        file<< runTime.timeName() << token::TAB << U << endl;
        runTime++;
    }

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
