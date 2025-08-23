scalar DiNum = -GREAT;

forAll(solidRegions, i)
{
    //- Note: do not use setRegionSolidFields.H to avoid double registering Cp
    //#include "setRegionSolidFields.H"
    const solidThermo& thermo = thermos[i];

    const volScalarField& cp = thermo.Cp();

    tmp<volScalarField> trho = thermo.rho();
    const volScalarField& rho = trho();

    DiNum = max
    (
        solidRegionDiffNo
        (
            solidRegions[i],
            runTime,
            rho*cp,
            thermo.isotropic() ? thermo.kappa() : mag(thermo.Kappa())()
        ),
        DiNum
    );

}
