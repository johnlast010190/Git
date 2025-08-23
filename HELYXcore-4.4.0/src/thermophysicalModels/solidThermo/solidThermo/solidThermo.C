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

#include "solidThermo/solidThermo.H"
#include "fvMesh/fvMesh.H"


/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(solidThermo, 0);
    defineRunTimeSelectionTable(solidThermo, objectRegistry);
    defineRunTimeSelectionTable(solidThermo, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidThermo::implementation::implementation
(
    const objectRegistry& obr,
    const word& phaseName
)
:
    p_(lookupOrConstruct(obr, "p", dimPressure)),
    rho_
    (
        IOobject
        (
            phasePropertyName("rho", phaseName),
            obr.time().timeName(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(obr),
        dimDensity
    )
{}


Foam::solidThermo::implementation::implementation
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName
)
:
    p_(lookupOrConstruct(obr, "p", dimPressure)),
    rho_
    (
        IOobject
        (
            phasePropertyName("rho", phaseName),
            obr.time().timeName(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(obr),
        dimDensity
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidThermo> Foam::solidThermo::New
(
    const objectRegistry& obr,
    const word& phaseName
)
{
    return basicThermo::New<solidThermo>(obr, phaseName);
}


Foam::autoPtr<Foam::solidThermo> Foam::solidThermo::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName
)
{
    return basicThermo::New<solidThermo>(obr, dict, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidThermo::~solidThermo()
{}


Foam::solidThermo::implementation::~implementation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::volScalarField& Foam::solidThermo::implementation::p()
{
    return p_;
}


const Foam::volScalarField& Foam::solidThermo::implementation::p() const
{
    return p_;
}

Foam::tmp<Foam::volScalarField> Foam::solidThermo::implementation::rho() const
{
    return rho_;
}


Foam::tmp<Foam::scalarField> Foam::solidThermo::implementation::rho
(
    const label patchi
) const
{
    return rho_.boundaryField()[patchi];
}


Foam::volScalarField& Foam::solidThermo::implementation::rho()
{
    return rho_;
}


// ************************************************************************* //
