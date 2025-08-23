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
    (c) 2019-2022 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "solidDisplacementThermo/solidDisplacementThermo.H"
#include "fvMesh/fvMesh.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(solidDisplacementThermo, 0);
}


Foam::materialTables& Foam::solidDisplacementThermo::matLookupOrConstruct
(
    const objectRegistry& obr,
    const dictionary& dict
)
{
    if (!obr.foundObject<objectRegistry>("materialModels"))
    {
        obr.store(new materialTables(obr, dict));
    }
    else if
    (
        !obr.subRegistry("materialModels").foundObject<materialTables>
        (
            "materialTables"
        )
    )
    {
        obr.store(new materialTables(obr, dict));
    }

    return
        obr.subRegistry
        (
            "materialModels"
        ).lookupObjectRef<materialTables>("materialTables");
}


void Foam::solidDisplacementThermo::correct()
{
    if (this->properties().lookup<word>("rhoModel") != "file")
    {
        rho_ = materials_(rhoModel::typeName)();
    }
    if (this->properties().lookup<word>("CpModel") != "file")
    {
        Cp_ = materials_(CpModel::typeName)();
    }
    if (this->properties().lookup<word>("kappaModel") != "file")
    {
        kappa_ = materials_(kappaModel::typeName)();
    }
    if (this->properties().lookup<word>("EModel") != "file")
    {
        E_ = materials_(EModel::typeName)();
    }
    if (this->properties().lookup<word>("nuModel") != "file")
    {
        nu_ = materials_(nuModel::typeName)();
    }
    if (this->properties().lookup<word>("alphavModel") != "file")
    {
        alphav_ = materials_(alphavModel::typeName)();
    }
}


void Foam::solidDisplacementThermo::loadModels(volScalarField& prop) const
{
    const word& propName = prop.name();
    if (this->properties().lookup<word>(propName + "Model") != "file")
    {
        materials_.addSpeciesAndSpeciesMixtures
        (
            propName + "Model",
            wordList({propName}),
            this->phaseName_,
            word::null,
            word::null,
            wordList({propName})
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidDisplacementThermo::solidDisplacementThermo
(
    const objectRegistry& obr,
    const word& phaseName
)
:
    solidThermo::composite(obr, phaseName),
    materials_(matLookupOrConstruct(obr, this->properties())),
    planeStress_(lookup("planeStress")),
    thermalStress_(lookup("thermalStress")),
    Cp_
    (
        IOobject
        (
            "Cp",
            obr.time().timeName(0),
            obr,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        this->mesh(obr),
        dimEnergy/dimMass/dimTemperature
    ),
    kappa_
    (
        IOobject
        (
            "kappa",
            obr.time().timeName(0),
            obr,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        this->mesh(obr),
        Cp_.dimensions()*dimensionSet(1, -1, -1, 0, 0)
    ),
    E_
    (
        IOobject
        (
            "E",
            obr.time().timeName(0),
            obr,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        this->mesh(obr),
        dimPressure
    ),
    nu_
    (
        IOobject
        (
            "nu",
            obr.time().timeName(0),
            obr,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        this->mesh(obr),
        dimless
    ),
    alphav_
    (
        IOobject
        (
            "alphav",
            obr.time().timeName(0),
            obr,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        this->mesh(obr),
        dimless/dimTemperature
    )
{
    loadModels(rho_);
    loadModels(Cp_);
    loadModels(kappa_);
    loadModels(E_);
    loadModels(nu_);
    loadModels(alphav_);

    // Update model pointers and dependency lists
    materials_.linkModelsAndUpdateTables();

    solidDisplacementThermo::correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidDisplacementThermo::~solidDisplacementThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::solidDisplacementThermo::rho() const
{
    return rho_;
}


Foam::tmp<Foam::scalarField> Foam::solidDisplacementThermo::rho
(
    const label patchi
) const
{
    return rho_.boundaryField()[patchi];
}


const Foam::volScalarField& Foam::solidDisplacementThermo::E() const
{
    return E_;
}


const Foam::scalarField& Foam::solidDisplacementThermo::E
(
    const label patchi
) const
{
    return E_.boundaryField()[patchi];
}


const Foam::volScalarField& Foam::solidDisplacementThermo::nu() const
{
    return nu_;
}


const Foam::scalarField& Foam::solidDisplacementThermo::nu
(
    const label patchi
) const
{
    return nu_.boundaryField()[patchi];
}


const Foam::volScalarField& Foam::solidDisplacementThermo::alphav() const
{
    return alphav_;
}


const Foam::scalarField& Foam::solidDisplacementThermo::alphav
(
    const label patchi
) const
{
    return alphav_.boundaryField()[patchi];
}


const Foam::volScalarField& Foam::solidDisplacementThermo::Cp() const
{
    return Cp_;
}


Foam::tmp<Foam::scalarField> Foam::solidDisplacementThermo::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    return Cp_.boundaryField()[patchi];
}


const Foam::volScalarField& Foam::solidDisplacementThermo::Cv() const
{
    return Cp_;
}


Foam::tmp<Foam::scalarField> Foam::solidDisplacementThermo::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    return Cp_.boundaryField()[patchi];
}


const Foam::volScalarField& Foam::solidDisplacementThermo::Cpv() const
{
    return Cp_;
}


Foam::tmp<Foam::scalarField> Foam::solidDisplacementThermo::Cpv
(
    const scalarField& T,
    const label patchi
) const
{
    return Cp_.boundaryField()[patchi];
}


const Foam::volScalarField& Foam::solidDisplacementThermo::kappa() const
{
    return kappa_;
}


Foam::tmp<Foam::scalarField> Foam::solidDisplacementThermo::kappa
(
    const label patchi
) const
{
    return kappa_.boundaryField()[patchi];
}


bool Foam::solidDisplacementThermo::read()
{
    return regIOobject::read();
}


// ************************************************************************* //