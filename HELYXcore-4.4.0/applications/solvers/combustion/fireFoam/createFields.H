Info<< "Reading thermophysical properties\n" << endl;

autoPtr<fluidMulticomponentThermo> pThermo(fluidMulticomponentThermo::New(mesh));
fluidMulticomponentThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

SLGThermo slgThermo(mesh, thermo);

speciesMassFractions& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

Info<< "Creating field rho\n" << endl;
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

volScalarField& p = thermo.p();

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

#include "cfdTools/general/include/readGravitationalAcceleration.H"
#include "cfdTools/general/include/readhRef.H"
#include "cfdTools/general/include/gh.H"

#include "cfdTools/compressible/compressibleCreatePhi.H"

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

Info<< "Creating combustion model\n" << endl;
autoPtr<combustionModel> combustion
(
    combustionModel::New(thermo, turbulence())
);


#include "cfdTools/general/include/readpRef.H"

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

mesh.schemes().setFluxRequired(p_rgh.name());

#include "phrghEqn.H"


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

IOdictionary additionalControlsDict
(
    IOobject
    (
        "additionalControls",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Switch solvePrimaryRegion
(
    additionalControlsDict.lookup("solvePrimaryRegion")
);

Switch solvePyrolysisRegion
(
    additionalControlsDict.lookupOrDefault<bool>("solvePyrolysisRegion", true)
);

#include "cfdTools/compressible/createDpdt.H"
#include "cfdTools/general/include/createK.H"
#include "cfdTools/general/include/createFvOptions.H"

Info<< "\nConstructing reacting cloud" << endl;
basicReactingCloud parcels("reactingCloud1", rho, U, g, slgThermo);

Info<< "\nConstructing surface film model" << endl;
autoPtr<regionModels::surfaceFilmModel> tsurfaceFilm
(
    regionModels::surfaceFilmModel::New(mesh, g)
);

Info<< "Creating pyrolysis model" << endl;
regionModels::pyrolysisModelCollection pyrolysis(mesh);

#include "solverObjects/radiation/createRadiationOption.H"