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
    (c) 2010-2016 Engys Ltd.
    (c) 2013-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "PhaseIncompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::PhaseIncompressibleTurbulenceModel<TransportModel>::
PhaseIncompressibleTurbulenceModel
(
    const word& type,
    const volScalarField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const TransportModel& transportModel,
    const word& propertiesName
)
:
    TurbulenceModel
    <
        volScalarField,
        geometricOneField,
        incompressibleTurbulenceModel,
        TransportModel
    >
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transportModel,
        propertiesName
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::autoPtr<Foam::PhaseIncompressibleTurbulenceModel<TransportModel>>
Foam::PhaseIncompressibleTurbulenceModel<TransportModel>::New
(
    const volScalarField& alpha,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const TransportModel& transportModel,
    const word& propertiesName
)
{
    autoPtr<PhaseIncompressibleTurbulenceModel> modelPtr
    (
        static_cast<PhaseIncompressibleTurbulenceModel*>
        (
            TurbulenceModel
            <
                volScalarField,
                geometricOneField,
                incompressibleTurbulenceModel,
                TransportModel
            >::New
            (
                alpha,
                geometricOneField(),
                U,
                alphaRhoPhi,
                phi,
                transportModel,
                propertiesName
            ).ptr()
        )
    );
    modelPtr->read();
    return modelPtr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TransportModel>
Foam::tmp<Foam::volScalarField>
Foam::PhaseIncompressibleTurbulenceModel<TransportModel>::pPrime() const
{
    return volScalarField::New
    (
        IOobject::groupName("pPrime", this->alphaRhoPhi_.group()),
        this->db_,
        this->mesh_,
        dimensionedScalar(dimPressure, 0)
    );
}


template<class TransportModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::PhaseIncompressibleTurbulenceModel<TransportModel>::pPrimef() const
{
    return surfaceScalarField::New
    (
        IOobject::groupName("pPrimef", this->alphaRhoPhi_.group()),
        this->db_,
        this->mesh_,
        dimensionedScalar(dimPressure, 0)
    );
}


template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::PhaseIncompressibleTurbulenceModel<TransportModel>::devReff() const
{
    return devRhoReff();
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::PhaseIncompressibleTurbulenceModel<TransportModel>::divDevReff
(
    volVectorField& U
) const
{
    return divDevRhoReff(U);
}

template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::PhaseIncompressibleTurbulenceModel<TransportModel>::divDevReff
(
    volVectorField& U,
    const word& scheme
) const
{
    return divDevRhoReff(U, scheme);
}


template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::PhaseIncompressibleTurbulenceModel<TransportModel>::
devRhoReff() const
{
    NotImplemented;
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::PhaseIncompressibleTurbulenceModel<TransportModel>::
divDevRhoReff
(
    volVectorField& U
) const
{
    NotImplemented;
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::PhaseIncompressibleTurbulenceModel<TransportModel>::
divDevRhoReff
(
    volVectorField& U,
    const word& scheme
) const
{
    NotImplemented;
}


// ************************************************************************* //
