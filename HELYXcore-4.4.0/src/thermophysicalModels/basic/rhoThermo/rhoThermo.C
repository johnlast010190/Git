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
    (c) 2011-2022 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "rhoThermo/rhoThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rhoThermo, 0);
    defineRunTimeSelectionTable(rhoThermo, objectRegistry);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rhoThermo::implementation::implementation
(
    const objectRegistry& obr,
    const word& phaseName
)
:
    rho_
    (
        IOobject
        (
            phasePropertyName("thermo:rho", phaseName),
            obr.time().timeName(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(obr),
        dimDensity
    )
{}


Foam::rhoThermo::implementation::implementation
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName
)
:
    rho_
    (
        IOobject
        (
            phasePropertyName("thermo:rho", phaseName),
            obr.time().timeName(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(obr),
        dimDensity
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::rhoThermo> Foam::rhoThermo::New
(
    const objectRegistry& obr,
    const word& phaseName
)
{
    return basicThermo::New<rhoThermo>(obr, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rhoThermo::~rhoThermo()
{}


Foam::rhoThermo::implementation::~implementation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::rhoThermo::implementation::rho() const
{
    return rho_;
}


Foam::tmp<Foam::scalarField> Foam::rhoThermo::implementation::rho
(
    const label patchi
) const
{
    return rho_.boundaryField()[patchi];
}


Foam::volScalarField& Foam::rhoThermo::implementation::rho()
{
    return rho_;
}


void Foam::rhoThermo::implementation::correctRho
(
    const Foam::volScalarField& deltaRho,
    const dimensionedScalar& rhoMin,
    const dimensionedScalar& rhoMax
)
{
    rho_ += deltaRho;
    rho_ = max(rho_, rhoMin);
    rho_ = min(rho_, rhoMax);
}


void Foam::rhoThermo::implementation::correctRho
(
    const volScalarField& deltaRho
)
{
    rho_ += deltaRho;
}


// ************************************************************************* //
