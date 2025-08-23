scalar nu2 = 0;

IOobject transportDictIO
(
    "transportProperties",
    runTime.constant(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::NO_WRITE,
    false
);

if (transportDictIO.typeHeaderOk<IOdictionary>())
{
    Info<< "Reading transportProperties\n" << endl;
    IOdictionary propertiesDict(transportDictIO);

    nu2 =
        dimensionedScalar
        (
            "nu",
            dimViscosity,
            propertiesDict.lookup("nu")
        ).value();
}
else
{
    Info<< "Reading materialProperties\n" << endl;
    IOdictionary materialProperties
    (
        IOobject
        (
            "materialProperties",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const scalar mu =
        materialProperties.subDict("muModelCoeffs")
        .lookup<scalar>("mu");

    const scalar rho =
        materialProperties.subDict("equationOfStateCoeffs")
        .lookup<scalar>("rho");

    nu2 = mu/rho;
}

Info<< "Reading viscoelastic properties\n" << endl;
IOdictionary turbulenceProperties
(
    IOobject
    (
        "turbulenceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
const dictionary& MaxwellCoeffs =
    turbulenceProperties.subDict("laminar").subDict("MaxwellCoeffs");
const scalar nu1 = MaxwellCoeffs.lookup<scalar>("nuM");
const scalar lambda = MaxwellCoeffs.lookup<scalar>("lambda");

const scalar rho = 1;

Info<< "Reading pressure gradient\n" << endl;
IOdictionary fvOptions
(
    IOobject
    (
        "fvOptions",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
const dictionary& gradPDict =
    fvOptions.subDict("momentumSource").subDict("injectionRateSuSp");
const scalar K =
    Tuple2<vector, scalar>(gradPDict.lookup("U")).first().x();

dictionary probes(IFstream(runTime.system()/"probes")());
const point location = pointField(probes.lookup("probeLocations"))[0];
