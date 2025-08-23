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
    (c) 2013-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "CompressibleTurbulenceModel/CompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::CompressibleTurbulenceModel<TransportModel>::
CompressibleTurbulenceModel
(
    const word& type,
    const geometricOneField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    TurbulenceModel
    <
        geometricOneField,
        volScalarField,
        compressibleTurbulenceModel,
        transportModel
    >
    (
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

template<class TransportModel>
Foam::autoPtr<Foam::CompressibleTurbulenceModel<TransportModel>>
Foam::CompressibleTurbulenceModel<TransportModel>::New
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
{
    autoPtr<CompressibleTurbulenceModel> pModel
    (
        static_cast<CompressibleTurbulenceModel*>(
        TurbulenceModel
        <
            geometricOneField,
            volScalarField,
            compressibleTurbulenceModel,
            transportModel
        >::New
        (
            geometricOneField(),
            rho,
            U,
            phi,
            phi,
            transport,
            propertiesName
        ).ptr())
    );
    pModel->read();
    return pModel;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class TransportModel>
Foam::tmp<Foam::fvBlockMatrix<Foam::vector>>
Foam::CompressibleTurbulenceModel<TransportModel>::bDivDevReff
(
    volVectorField& U
) const
{
    return bDivDevRhoReff(U);
}

template<class TransportModel>
Foam::tmp<Foam::fvBlockMatrix<Foam::vector>>
Foam::CompressibleTurbulenceModel<TransportModel>::
bDivDevRhoReff
(
    volVectorField& U
) const
{
    NotImplemented;
}


template<class TransportModel>
Foam::tmp<Foam::fvBlockMatrix<Foam::vector>>
Foam::CompressibleTurbulenceModel<TransportModel>::bDivDevReff
(
    volVectorField& U,
    const word& scheme
) const
{
    return bDivDevRhoReff(U, scheme);
}

template<class TransportModel>
Foam::tmp<Foam::fvBlockMatrix<Foam::vector>>
Foam::CompressibleTurbulenceModel<TransportModel>::
bDivDevRhoReff
(
    volVectorField& U,
    const word& scheme
) const
{
    NotImplemented;
}



// ************************************************************************* //
