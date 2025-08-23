/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
|    o     o     |  ENGYS Ltd. <http://engys.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of HELYXcore.
    HELYXcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    HELYXcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HELYXcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with HELYXcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2021-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "HerschelBulkley.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "finiteVolume/fvc/fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixtureViscosityModels
{
    defineTypeNameAndDebug(HerschelBulkley, 0);

    addToRunTimeSelectionTable
    (
        mixtureViscosityModel,
        HerschelBulkley,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixtureViscosityModels::HerschelBulkley::HerschelBulkley
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word modelName
)
:
    mixtureViscosityModel(name, viscosityProperties, U, phi),
    HerschelBulkleyCoeffs_
    (
        viscosityProperties.optionalSubDict(modelName + "Coeffs")
    ),
    k_
    (
        "coeff",
        dimensionSet(1, -1, -1, 0, 0),
        HerschelBulkleyCoeffs_.lookupOrDefault<scalar>("coeff", 0.1)
    ),
    n_
    (
        "exponent",
        dimless,
        HerschelBulkleyCoeffs_.lookupOrDefault<scalar>("exponent", 0.6)
    ),
    tau0_
    (
        "tau0",
        dimensionSet(1, -1, -1, 0, 0),
        HerschelBulkleyCoeffs_.lookupOrDefault<scalar>("tau0", 0.0)
    ),
    mu0_
    (
        "mu0",
        dimensionSet(1, -1, -1, 0, 0),
        HerschelBulkleyCoeffs_.lookup("mu0")
    ),
    muMax_
    (
        "muMax",
        dimensionSet(1, -1, -1, 0, 0),
        HerschelBulkleyCoeffs_.lookup("muMax")
    ),
    fixedParam_
    (
        HerschelBulkleyCoeffs_.lookupOrDefault<bool>("fixedParam", false)
    ),
    isDrySolids_
    (
        HerschelBulkleyCoeffs_.lookupOrDefault<bool>("drySolidsInput", true)
    ),
    rhod_
    (
        "rhod",
        dimless,
        HerschelBulkleyCoeffs_.lookup("rhod")
    ),
    alpha_
    (
        U.mesh().lookupObject<volScalarField>
        (
            HerschelBulkleyCoeffs_.lookupOrDefault<word>("alpha", "alphad")
        )
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::mixtureViscosityModels::HerschelBulkley::mu(const volScalarField& muc) const
{
    dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);

    tmp<volScalarField> sr(strainRate());

    if (U_.mesh().time().outputTime() && debug)
    {
        volScalarField strainRate
        (
            IOobject
            (
                "strainRate",
                U_.time().timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sr()
        );
        strainRate.write();
    }

    return
    (
        min
        (
            muMax_,
            max
            (
                mu0_,
                (tau0() + k()*pow(tone*sr(), n())) * rtone
               /(max(sr(), dimensionedScalar(dimless/dimTime, VSMALL)))
            )
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::mixtureViscosityModels::HerschelBulkley::strainRate() const
{
    return sqrt(2.0)*mag(symm(fvc::grad(U_)));
}


Foam::tmp<Foam::volScalarField>
Foam::mixtureViscosityModels::HerschelBulkley::n() const
{
    tmp<volScalarField> exponent(volScalarField::New("n", U_.mesh(), n_));
    if (!fixedParam_)
    {
        exponent.ref().forceAssign(min(0.6019*pow(X(), -0.172), scalar(1.0)));
    }
    return exponent;
}


Foam::tmp<Foam::volScalarField>
Foam::mixtureViscosityModels::HerschelBulkley::k() const
{
    tmp<volScalarField> coeff(volScalarField::New("k", U_.mesh(), k_));
    if (!fixedParam_)
    {
        dimensionedScalar a("a", k_.dimensions(), scalar(0.0934));
        coeff.ref().forceAssign(a*pow(X(), 1.681));
    }
    return coeff;
}


Foam::tmp<Foam::volScalarField>
Foam::mixtureViscosityModels::HerschelBulkley::tau0() const
{
    tmp<volScalarField> tau0(volScalarField::New("tau0", U_.mesh(), tau0_));
    if (!fixedParam_)
    {
        dimensionedScalar a("a", tau0_.dimensions(), scalar(0.021));
        tau0.ref().forceAssign(a*pow(X(), 2.6891));
    }
    return tau0;
}


Foam::tmp<Foam::volScalarField>
Foam::mixtureViscosityModels::HerschelBulkley::X() const
{
    tmp<volScalarField> alpha(max(alpha_, VSMALL));

    if (isDrySolids_)
    {
        return alpha()*100.0;
    }

    return alpha*rhod_;
}


bool Foam::mixtureViscosityModels::HerschelBulkley::read
(
    const dictionary& viscosityProperties
)
{
    mixtureViscosityModel::read(viscosityProperties);

    HerschelBulkleyCoeffs_ =
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    return true;
}


// ************************************************************************* //
