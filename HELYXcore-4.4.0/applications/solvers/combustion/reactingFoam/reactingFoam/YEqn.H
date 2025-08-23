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
              - fvm::laplacian(turbulence->muEff(), Yi)
             ==
                combustion->R(Yi)
              + fvOptions(rho, Yi)
            );

            // Adjusted order of addition to match USF
            YiEqn += mvConvection->fvmDiv(phi, Yi);
            YiEqn.relax("Yi");
            fvOptions.constrain(YiEqn);
            YiEqn.solve(mesh.solution().solver("Yi"));
            fvOptions.correct(Yi);
        }
    }
    composition.normalise();
}
