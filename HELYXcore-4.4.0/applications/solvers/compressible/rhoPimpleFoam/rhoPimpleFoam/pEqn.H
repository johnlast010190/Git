if (!pimple.SIMPLErho())
{
    rho = thermo.rho();
}

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField p0("p0", p);
tmp<surfaceScalarField> rhof(fvc::interpolate(rho));

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", rhof()*fvc::interpolate(rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

surfaceScalarField phiHbyA
(
    "phiHbyA",
    rhof()*fvc::flux(HbyA)
);
if (ddtPhiCorr)
{
    phiHbyA += rhorAUf*fvc::ddtCorr(rho, U, phi);
}

fvOptions.makeRelative(rhof(), phiHbyA);


// Update the pressure BCs to ensure flux consistency
constrainPressure(p, rho, U, phiHbyA, rhorAUf, fvOptions);

if (pimple.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        (fvc::interpolate(psi)/rhof())*phiHbyA
    );

    phiHbyA -= fvc::interpolate(psi*p)*phiHbyA/rhof();

    fvScalarMatrix pDDtEqn
    (
        fvc::ddt(rho) + psi*correction(fvm::ddt(p))
      + fvc::div(phiHbyA) + fvm::div(phid, p)
     ==
        fvOptions(psi, p, rho.name())
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAUf, p));

        // Relax the pressure equation to ensure diagonal-dominance
        pEqn.relax();

        pEqn.solve(mesh.solution().solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }
}
else
{
    fvScalarMatrix pDDtEqn
    (
      - fvOptions(psi, p, rho.name())
      + fvc::div(phiHbyA)
    );

    // Different addition order to match USF
    pDDtEqn += fvc::ddt(rho);
    pDDtEqn += psi*correction(fvm::ddt(p));

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAUf, p));

        pEqn.solve(mesh.solution().solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }
}

#include "cfdTools/compressible/rhoEqn.H"
#include "cfdTools/compressible/compressibleContinuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);
K = 0.5*magSqr(U);

if (pressureControl.limit(p))
{
    p.correctBoundaryConditions();
}

thermo.correctRho(psi*(p - p0), rhoMin, rhoMax);
rho = thermo.rho();

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}
