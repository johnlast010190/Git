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
    (c) 2010-2012 Engys Ltd.
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "viscosityModels/powerLaw/powerLaw.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "finiteVolume/fvc/fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(powerLaw, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        powerLaw,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::powerLaw::calcNu() const
{
    return max
    (
        nuMin_,
        min
        (
            nuMax_,
            k_*pow
            (
                max
                (
                    dimensionedScalar(dimTime, 1)*strainRate(),
                    dimensionedScalar(dimless, VSMALL)
                ),
                n_.value() - scalar(1.0)
            )
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::powerLaw::powerLaw
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    powerLawCoeffs_(viscosityProperties.optionalSubDict(typeName + "Coeffs")),
    k_("k", dimViscosity, powerLawCoeffs_),
    n_("n", dimless, powerLawCoeffs_),
    nuMin_("nuMin", dimViscosity, powerLawCoeffs_),
    nuMax_("nuMax", dimViscosity, powerLawCoeffs_),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::viscosityModels::powerLaw::strainRate() const
{
    return (sqrt(2.0)*mag(symm((fvc::grad(U_)))));
}

bool Foam::viscosityModels::powerLaw::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    powerLawCoeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");

    powerLawCoeffs_.lookup("k") >> k_;
    powerLawCoeffs_.lookup("n") >> n_;
    powerLawCoeffs_.lookup("nuMin") >> nuMin_;
    powerLawCoeffs_.lookup("nuMax") >> nuMax_;

    return true;
}


// ************************************************************************* //
