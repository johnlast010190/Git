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
    (c) 2012-2022 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fluidThermo/fluidThermo.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(fluidThermo, 0);
    defineRunTimeSelectionTable(fluidThermo, objectRegistry);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidThermo::implementation::implementation
(
    const objectRegistry& obr,
    const word& phaseName
)
:
    p_(lookupOrConstruct(obr, "p", dimPressure)),
    psi_
    (
        IOobject
        (
            phasePropertyName("psi", phaseName),
            obr.time().timeName(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(obr),
        dimensionSet(0, -2, 2, 0, 0)
    ),
    mu_
    (
        IOobject
        (
            phasePropertyName("mu", phaseName),
            obr.time().timeName(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(obr),
        dimensionSet(1, -1, -1, 0, 0)
    )
{}


Foam::fluidThermo::implementation::implementation
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName
)
:
    p_(lookupOrConstruct(obr, "p", dimPressure)),
    psi_
    (
        IOobject
        (
            phasePropertyName("psi", phaseName),
            obr.time().timeName(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(obr),
        dimensionSet(0, -2, 2, 0, 0)
    ),
    mu_
    (
        IOobject
        (
            phasePropertyName("mu", phaseName),
            obr.time().timeName(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(obr),
        dimensionSet(1, -1, -1, 0, 0)
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fluidThermo> Foam::fluidThermo::New
(
    const objectRegistry& obr,
    const word& phaseName
)
{
    return basicThermo::New<fluidThermo>(obr, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidThermo::~fluidThermo()
{}


Foam::fluidThermo::implementation::~implementation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::fluidThermo::nu() const
{
    return mu()/rho();
}


Foam::tmp<Foam::scalarField> Foam::fluidThermo::nu(const label patchi) const
{
    return mu().boundaryField()[patchi]/rho(patchi);
}


Foam::volScalarField& Foam::fluidThermo::implementation::p()
{
    return p_;
}


const Foam::volScalarField& Foam::fluidThermo::implementation::p() const
{
    return p_;
}


const Foam::volScalarField& Foam::fluidThermo::implementation::psi() const
{
    return psi_;
}


const Foam::volScalarField& Foam::fluidThermo::implementation::mu() const
{
    return mu_;
}


// ************************************************************************* //
