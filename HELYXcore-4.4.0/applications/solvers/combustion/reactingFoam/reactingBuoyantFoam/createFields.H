#include "cfdTools/general/include/createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;
autoPtr<fluidMulticomponentThermo> pThermo(fluidMulticomponentThermo::New(mesh));
fluidMulticomponentThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

speciesMassFractions& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
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


volScalarField& p = thermo.p();

#include "cfdTools/general/include/readGravitationalAcceleration.H"
#include "cfdTools/general/include/readhRef.H"
#include "cfdTools/general/include/gh.H"

#include "cfdTools/compressible/compressibleCreatePhi.H"

mesh.schemes().setFluxRequired(p.name());

Info<< "Creating turbulence model.\n" << nl;
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

Info<< "Creating reaction model\n" << endl;
autoPtr<combustionModel> combustion
(
    combustionModel::New(thermo, turbulence())
);


Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Force p_rgh to be consistent with p
p_rgh = p - rho*gh;

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

#include "cfdTools/compressible/createDpdt.H"

#include "cfdTools/general/include/createK.H"

#include "cfdTools/general/include/createFvOptions.H"
