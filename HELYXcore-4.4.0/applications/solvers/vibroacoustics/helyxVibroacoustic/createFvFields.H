Info<< "\nCreating Fv Fields" << endl;

volScalarField wV
(
    IOobject
    (
        "wV",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimLength, 0)
);

volScalarField loadV
(
    IOobject
    (
        "loadV",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimensionSet(1, -1, -2, 0, 0, 0, 0), 0)
);

volScalarField pInt
(
    IOobject
    (
        "pInt",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
