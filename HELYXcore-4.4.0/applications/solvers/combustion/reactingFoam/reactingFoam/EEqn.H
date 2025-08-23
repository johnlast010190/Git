{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        mvConvection->fvmDiv(phi, he)
      + fvc::div(phi, K)
     ==
        combustion->Qdot()
      + fvOptions(rho, he)
    );

    // Match to USF
    EEqn += fvm::ddt(rho, he) + fvc::ddt(rho, K);
    EEqn +=
        (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        );
    EEqn -= fvm::laplacian(turbulence->alphaEff(), he);

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
