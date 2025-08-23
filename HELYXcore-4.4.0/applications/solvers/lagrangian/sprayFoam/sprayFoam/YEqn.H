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
              - fvm::laplacian(turbulence->muEff(), Yi)
             ==
                parcels.SYi(i, Yi)
              + combustion->R(Yi)
              + fvOptions(rho, Yi)
            );

            YEqn.relax();
            fvOptions.constrain(YEqn);
            YEqn.solve(mesh.solution().solver("Yi"));
            fvOptions.correct(Yi);
        }
    }
    composition.normalise();
}
