tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.schemes().divScheme("div(phi,Yi_h)")
    )
);
{
    combustion->correct();
    forAll(Y, i)
    {
        if (composition.solve(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turbulence->alphaEff(), Yi)
              ==
                parcels.SYi(i, Yi)
              + surfaceFilm.Srho(i)
              + combustion->R(Yi)
              + fvOptions(rho, Yi)
            );

            YiEqn.relax();
            fvOptions.constrain(YiEqn);
            YiEqn.solve(mesh.solution().solver("Yi"));
            fvOptions.correct(Yi);
        }
    }
    composition.normalise();

    radiation.correct();

    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(turbulence->alphaEff(), he)
     ==
        combustion->Qdot()
      + parcels.Sh(he)
      + surfaceFilm.Sh()
      + fvOptions(rho, he)
    );

    EEqn.relax();
    fvOptions.constrain(EEqn);
    EEqn.solve();
    fvOptions.correct(he);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
