    if (mesh.nCells() != 1)
    {
        FatalErrorInFunction
            << "Solver only applicable to single cell cases"
            << exit(FatalError);
    }

    Info<< "Reading initial conditions.\n" << endl;
    IOdictionary initialConditions
    (
        IOobject
        (
            "initialConditions",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    scalar p0 = initialConditions.lookup<scalar>("p");
    scalar T0 = initialConditions.lookup<scalar>("T");

    // write base thermo fields - not registered since will be re-read by
    // thermo package

    Info<< "Creating base fields for time " << runTime.timeName() << endl;
    {
        volScalarField Ydefault
        (
            IOobject
            (
                "Ydefault",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar(dimless, 1)
        );

        Ydefault.write();

        volScalarField p
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar(dimPressure, p0)
        );

        p.write();

        volScalarField T
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar(dimTemperature, T0)
        );

        T.write();
    }

    Info<< "Reading thermophysical properties\n" << endl;
    autoPtr<fluidMulticomponentThermo> pThermo(fluidMulticomponentThermo::New(mesh));
    fluidMulticomponentThermo& thermo = pThermo();
    thermo.validate(args.executable(), "h");

    autoPtr<basicChemistryModel> pChemistry(basicChemistryModel::New(thermo));

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

    volScalarField Rspecific
    (
        IOobject
        (
            "Rspecific",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, 0)
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector(dimVelocity, Zero)
    );

    #include "cfdTools/incompressible/createPhi.H"

    OFstream post(args.path()/"chemFoam.out");
    post<< "# Time" << token::TAB << "Temperature [K]" << token::TAB
        << "Pressure [Pa]" << endl;
