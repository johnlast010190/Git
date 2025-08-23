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
    (c) 2013-2019 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "TurbulenceModel.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    class Alpha,
    class Rho,
    class BasicTurbulenceModel,
    class TransportModel
>
Foam::TurbulenceModel<Alpha, Rho, BasicTurbulenceModel, TransportModel>::
TurbulenceModel
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicTurbulenceModel
    (
        rho,
        U,
        alphaRhoPhi,
        phi,
        propertiesName
    ),
    alpha_(alpha),
    transport_(transport)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template
<
    class Alpha,
    class Rho,
    class BasicTurbulenceModel,
    class TransportModel
>
Foam::autoPtr
<
    Foam::TurbulenceModel<Alpha, Rho, BasicTurbulenceModel, TransportModel>
>
Foam::TurbulenceModel<Alpha, Rho, BasicTurbulenceModel, TransportModel>::New
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
{
    // get model name, but do not register the dictionary
    // otherwise it is registered in the database twice
    const word modelType
    (
        IOdictionary
        (
            IOobject
            (
                IOobject::groupName(propertiesName, alphaRhoPhi.group()),
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("simulationType")
    );

    Info<< "Selecting turbulence model type " << modelType << endl;

    const auto ctor =
        ctorTableLookup
        (
            "TurbulenceModel type",
            dictionaryConstructorTable_(),
            modelType
        );
    return autoPtr<TurbulenceModel>
    (
        ctor(alpha, rho, U, alphaRhoPhi, phi, transport, propertiesName)
    );
}


// ************************************************************************* //
