{
    if (!pimple.SIMPLErho())
    {
        rho = thermo.rho();
    }

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution
    const volScalarField psip0(psi*p);

    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            fvc::flux(rho*HbyA)
          + rhorAUf*fvc::ddtCorr(rho, U, phi)
        )
      + phig
    );

    {
        tmp<surfaceScalarField> rhof(fvc::interpolate(rho));
        fvOptions.makeRelative(rhof, phiHbyA);
    }

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, fvOptions);

    fvScalarMatrix p_rghDDtEqn
    (
        fvc::ddt(rho) + psi*correction(fvm::ddt(p_rgh))
      + fvc::div(phiHbyA)
     ==
        parcels.Srho()
      + surfaceFilm.Srho()
      + fvOptions(psi, p_rgh, rho.name())
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            p_rghDDtEqn
          - fvm::laplacian(rhorAUf, p_rgh)
        );

        p_rghEqn.solve(mesh.solution().solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + p_rghEqn.flux();
            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            U = HbyA + rAU*fvc::reconstruct((p_rghEqn.flux() + phig)/rhorAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
            K = 0.5*magSqr(U);
        }
    }

    p = p_rgh + rho*gh;

    bool limitedp = pressureControl.limit(p);

    if (limitedp)
    {
        p_rgh = p - rho*gh;
    }

    // Thermodynamic density update
    thermo.correctRho(psi*p - psip0);

    if (limitedp)
    {
        rho = thermo.rho();
    }

    #include "rhoEqn.H" // NOTE: flux and time scales now inconsistent
    #include "cfdTools/compressible/compressibleContinuityErrs.H"

    if (pimple.SIMPLErho())
    {
        rho = thermo.rho();
    }

    if (thermo.dpdt())
    {
        dpdt = fvc::ddt(p);
    }
}
