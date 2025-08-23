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

#include "phasePressureModel.H"
#include "twoPhaseSystem.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RASModels::phasePressureModel::phasePressureModel
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& phase,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity
    <
        RASModel<EddyDiffusivity<phaseCompressibleTurbulenceModel>>
    >
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        phase,
        propertiesName
    ),
    alphaMax_(coeffDict_.lookup<scalar>("alphaMax")),
    preAlphaExp_(coeffDict_.lookup<scalar>("preAlphaExp")),
    expMax_(coeffDict_.lookup<scalar>("expMax")),
    g0_
    (
        "g0",
        dimensionSet(1, -1, -2, 0, 0),
        coeffDict_.lookup("g0")
    )
{
    nut_.forceAssign(dimensionedScalar(nut_.dimensions(), 0));

    if (type == typeName)
    {
        printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RASModels::phasePressureModel::~phasePressureModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::RASModels::phasePressureModel::read()
{
    if
    (
        eddyViscosity
        <
            RASModel<EddyDiffusivity<phaseCompressibleTurbulenceModel>>
        >::read()
    )
    {
        alphaMax_ = coeffDict().lookup<scalar>("alphaMax");
        preAlphaExp_ = coeffDict().lookup<scalar>("preAlphaExp");
        expMax_ = coeffDict().lookup<scalar>("expMax");
        g0_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::phasePressureModel::k() const
{
    NotImplemented;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::phasePressureModel::epsilon() const
{
    NotImplemented;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::phasePressureModel::omega() const
{
    NotImplemented;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::phasePressureModel::R() const
{
    return volSymmTensorField::New
    (
        IOobject::groupName("R", U_.group()),
        mesh_,
        dimensioned<symmTensor>
        (
            "R",
            dimensionSet(0, 2, -2, 0, 0),
            Zero
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::phasePressureModel::pPrime() const
{
    tmp<volScalarField> tpPrime
    (
        g0_
       *min
        (
            exp(preAlphaExp_*(alpha_ - alphaMax_)),
            expMax_
        )
    );

    volScalarField::Boundary& bpPrime =
        tpPrime.ref().boundaryFieldRef();

    forAll(bpPrime, patchi)
    {
        if (!bpPrime[patchi].coupled())
        {
            bpPrime[patchi].forceAssign(0);
        }
    }

    return tpPrime;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::RASModels::phasePressureModel::pPrimef() const
{
    tmp<surfaceScalarField> tpPrime
    (
        g0_
       *min
        (
            exp(preAlphaExp_*(fvc::interpolate(alpha_) - alphaMax_)),
            expMax_
        )
    );

   surfaceScalarField::Boundary& bpPrime =
       tpPrime.ref().boundaryFieldRef();

    forAll(bpPrime, patchi)
    {
        if (!bpPrime[patchi].coupled())
        {
            bpPrime[patchi].forceAssign(0);
        }
    }

    return tpPrime;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::phasePressureModel::devRhoReff() const
{
    return volSymmTensorField::New
    (
        IOobject::groupName("devRhoReff", U_.group()),
        mesh_,
        dimensioned<symmTensor>
        (
            "R",
            rho_.dimensions()*dimensionSet(0, 2, -2, 0, 0),
            Zero
        )
    );
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::RASModels::phasePressureModel::divDevRhoReff
(
    volVectorField& U
) const
{
    return tmp<fvVectorMatrix>
    (
        new fvVectorMatrix
        (
            U,
            rho_.dimensions()*dimensionSet(0, 4, -2, 0, 0)
        )
    );
}


void Foam::RASModels::phasePressureModel::correct()
{}


// ************************************************************************* //
