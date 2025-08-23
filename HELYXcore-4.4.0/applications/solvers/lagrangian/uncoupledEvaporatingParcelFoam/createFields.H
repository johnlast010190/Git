/*---------------------------------------------------------------------------*\
| Modified 2010-2016 Copyright (C) Engys Ltd                                  |
\*---------------------------------------------------------------------------*/

    Info<< "Reading thermophysical properties\n" << endl;
    autoPtr<fluidMulticomponentThermo> pThermo(fluidMulticomponentThermo::New(mesh));
    fluidMulticomponentThermo& thermo = pThermo();
    thermo.validate(args.executable(), "h", "e");

    SLGThermo slgThermo(mesh, thermo);

    speciesMassFractions& composition = thermo.composition();
    PtrList<volScalarField>& Y = composition.Y();

    volScalarField& p = thermo.p();

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );

    Info<< "\nReading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "cfdTools/compressible/compressibleCreatePhi.H"

    dimensionedScalar rhoMax
    (
        "rhoMax",
        dimDensity,
        mesh.solution().dict().subDict("PIMPLE")
       .lookupOrDefault<scalar>("rhoMax", GREAT)
    );

    dimensionedScalar rhoMin
    (
        "rhoMin",
        dimDensity,
        mesh.solution().dict().subDict("PIMPLE")
       .lookupOrDefault<scalar>("rhoMin", -GREAT)
    );

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    Info<< "Creating field dpdt\n" << endl;
    volScalarField dpdt
    (
        IOobject
        (
            "dpdt",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(p.dimensions()/dimTime, 0)
    );

    Info<< "Creating field kinetic energy K\n" << endl;
    volScalarField K("K", 0.5*magSqr(U));

    Info<< "Creating multi-variate interpolation scheme\n" << endl;
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

    forAll(Y, i)
    {
        fields.add(Y[i]);
    }
    fields.add(thermo.he());

    volScalarField dQ
    (
        IOobject
        (
            "dQ",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimTime, 0)
    );

    volVectorField E
    (
        IOobject
        (
            "E",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
