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
    (c) 2013-2014 OpenFOAM Foundation
    (c) 2017-2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "incompressibleTurbulenceModel.H"
#include "turbulentTransportModels/derivedFvPatchFields/wallFunctions/alphatWallFunctions/alphatKaderWallFunction/alphatKaderWallFunctionFvPatchScalarField.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleTurbulenceModel, 0);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::incompressibleTurbulenceModel::alphatReadIfPresent
(
    const word& fieldName,
    const fvMesh& mesh
) const
{
    IOobject alphatHeader
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (alphatHeader.typeHeaderOk<volScalarField>(true))
    {
        tmp<volScalarField> talphat(new volScalarField(alphatHeader, mesh));

        return talphat;
    }
    else
    {
        return tmp<volScalarField>(nullptr);
    }
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleTurbulenceModel::alphatAutoCreate
(
    const word& fieldName,
    const fvMesh& mesh
) const
{
    Info<< "--> Creating " << fieldName
        << " to employ run-time selectable wall functions"
        << endl;

    const fvBoundaryMesh& bm = mesh.boundary();

    wordList alphatBoundaryTypes(bm.size());

    forAll(bm, patchi)
    {
        if (isA<wallFvPatch>(bm[patchi]))
        {
            alphatBoundaryTypes[patchi] =
                incompressible::
                alphatKaderWallFunctionFvPatchScalarField::typeName;
        }
        else
        {
            alphatBoundaryTypes[patchi] =
                calculatedFvPatchField<scalar>::typeName;
        }
    }

    tmp<volScalarField> alphat
    (
        new volScalarField
        (
            IOobject
            (
                fieldName,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimArea/dimTime, 0),
            alphatBoundaryTypes
        )
    );

    return alphat;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleTurbulenceModel::incompressibleTurbulenceModel
(
    const geometricOneField&,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const word& propertiesName
)
:
    turbulenceModel(U, alphaRhoPhi, phi, propertiesName),
    alphat_(nullptr)
{
    const word alphatName(IOobject::groupName("alphat", alphaRhoPhi.group()));
    tmp<volScalarField> alphaTemp(alphatReadIfPresent(alphatName, mesh_));

    if (alphaTemp.valid())
    {
        alphat_.set(alphaTemp.ptr());
    }
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::incompressibleTurbulenceModel::~incompressibleTurbulenceModel()
{}


const Foam::volScalarField& Foam::incompressibleTurbulenceModel::alphat() const
{
    // Create alphat when needed
    if (!alphat_.valid())
    {
        const word alphatName(IOobject::groupName("alphat", this->group()));
        alphat_.set(alphatAutoCreate(alphatName, mesh_).ptr());
        alphat_() = nut()/Prt();
        alphat_().correctBoundaryConditions();
    }

    // Backward compatibility
    if (alphat_().dimensions() == dimMass/dimLength/dimTime)
    {
        alphat_() /= rho();
    }

    return alphat_();
}


const Foam::tmp<Foam::volScalarField>
Foam::incompressibleTurbulenceModel::alphaEff() const
{
    return (alphat() + alphaLam());
}


const Foam::volScalarField&
Foam::incompressibleTurbulenceModel::mu() const
{
    nu_.reset(nu().ptr());
    return nu_();
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleTurbulenceModel::mut() const
{
    return nut();
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleTurbulenceModel::mut(const label patchi) const
{
    return nut(patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleTurbulenceModel::muEff() const
{
    return nuEff();
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleTurbulenceModel::muEff(const label patchi) const
{
    return nuEff(patchi);
}


void Foam::incompressibleTurbulenceModel::correct()
{
    turbulenceModel::correct();

    // this will be superfluous for laminar model, alphat=0
    if (alphat_.valid())
    {
        alphat_() = nut()/Prt();
        alphat_().correctBoundaryConditions();
    }
}


Foam::tmp<Foam::scalarField> Foam::incompressibleTurbulenceModel::uWallCoeffs
(
    const label& patchi
)
{
    return muEff(patchi)/this->y()[patchi];
}

// ************************************************************************* //
