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
    (c) 2015-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "ThermalDiffusivity.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::ThermalDiffusivity<BasicTurbulenceModel>::ThermalDiffusivity
(
    const word& type,
    const alphaField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicTurbulenceModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::autoPtr<Foam::ThermalDiffusivity<BasicTurbulenceModel>>
Foam::ThermalDiffusivity<BasicTurbulenceModel>::New
(
    const alphaField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
{
    return autoPtr<ThermalDiffusivity>
    (
        static_cast<ThermalDiffusivity*>
        (
            BasicTurbulenceModel::New
            (
                alpha,
                rho,
                U,
                alphaRhoPhi,
                phi,
                transport,
                propertiesName
            ).ptr()
        )
    );
}


template<class BasicTurbulenceModel>
Foam::autoPtr<Foam::ThermalDiffusivity<BasicTurbulenceModel>>
Foam::ThermalDiffusivity<BasicTurbulenceModel>::New
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
{
    return autoPtr<ThermalDiffusivity>
    (
        static_cast<ThermalDiffusivity*>
        (
            BasicTurbulenceModel::New
            (
                rho,
                U,
                phi,
                transport,
                propertiesName
            ).ptr()
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::ThermalDiffusivity<BasicTurbulenceModel>::alphat() const
{
    return volScalarField::New
    (
        IOobject::groupName("alphat", this->alphaRhoPhi_.group()),
        this->db_,
        this->mesh_,
        dimensionedScalar(dimDensity*dimViscosity, 0)
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::scalarField>
Foam::ThermalDiffusivity<BasicTurbulenceModel>::alphat
(
    const label patchi
) const
{
    return tmp<scalarField>
    (
        new scalarField(this->mesh_.boundary()[patchi].size(), 0.0)
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::ThermalDiffusivity<BasicTurbulenceModel>::Dmt() const
{
    return volScalarField::New
    (
        IOobject::groupName("Dmt", this->alphaRhoPhi_.group()),
        this->db_,
        this->mesh_,
        dimensionedScalar(dimDensity*dimViscosity, 0)
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::scalarField>
Foam::ThermalDiffusivity<BasicTurbulenceModel>::Dmt
(
    const label patchi
) const
{
    return tmp<scalarField>
    (
        new scalarField(this->mesh_.boundary()[patchi].size(), 0.0)
    );
}


// ************************************************************************* //
