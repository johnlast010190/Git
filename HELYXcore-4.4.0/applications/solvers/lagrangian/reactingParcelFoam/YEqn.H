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

            fvScalarMatrix YEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turbulence->alphaEff(), Yi)
             ==
                parcels.SYi(i, Yi)
              + fvOptions(rho, Yi)
              + combustion->R(Yi)
              + surfaceFilm.Srho(i)
            );

            YEqn.relax();
            fvOptions.constrain(YEqn);
            YEqn.solve(mesh.solution().solver("Yi"));
            fvOptions.correct(Yi);
        }
    }
    composition.normalise();
}
