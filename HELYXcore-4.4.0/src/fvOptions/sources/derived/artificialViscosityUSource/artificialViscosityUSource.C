/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX® : Professional Open-source CFD
|   o   O   o    |
|    o     o     |  ENGYS Ltd. <http://engys.com/>
|       o        |
\*---------------------------------------------------------------------------

License
    This file is part of HELYXcore.
    HELYXcore is based on OpenFOAM® <http://www.openfoam.org/>.

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
    (c) 2018-2024 Engys Ltd.

Author
    2018. Nikolaos Magoulas (Engys Ltd.). All rights reserved.

\*---------------------------------------------------------------------------*/

#include "artificialViscosityUSource.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "turbulenceModel.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchField.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(artificialViscosityUSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        artificialViscosityUSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::artificialViscosityUSource::artificialViscosityUSource
(
    const word& name,
    const word& modelType,
    const dictionary& optionDict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, optionDict, obr),
    sensor_
    (
        sensor<vector>::New(mesh_, optionDict.subDict("sensor"))
    ),
    viscosityCoeff_
    (
        optionDict.lookupOrDefault<scalar>("viscosityCoeff", 0.0)
    ),
    writeViscosity_
    (
        optionDict.lookupOrDefault<Switch>("writeViscosity", false)
    ),
    maxViscosity_
    (
        optionDict.lookupOrDefault<scalar>("maxViscosity", 1.0/SMALL)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::artificialViscosityUSource::sourceFields(wordList& fieldNames)
{
    fieldNames = coeffs_.lookup<wordList>("fields");
}


void Foam::fv::artificialViscosityUSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    dimensionedScalar viscCoeff
    (
        "visc",
        eqn.psi().dimensions()*dimLength*dimLength/dimTime,
        viscosityCoeff_
    );

    volScalarField viscosity
    (
        IOobject
        (
            "artificialViscosity_" + name(),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        viscCoeff,
        zeroGradientFvPatchField<scalar>::typeName
    );
    viscosity.primitiveFieldRef() *= sensor_().valueField();
    viscosity.primitiveFieldRef() =
        min(viscosity.primitiveField(), maxViscosity_);
    viscosity.correctBoundaryConditions();

    if (mesh_.time().outputTime())
    {
        sensor_().writeField();

        if (writeViscosity_ == true)
        {
            viscosity.write();
        }
    }

    eqn +=
        fvm::laplacian(viscosity, eqn.psi())
      + fvc::div
        (
            viscosity*dev2(T(fvc::grad(eqn.psi()))),
            "div((nu*dev2(T(grad(U)))))"
        );
}

void Foam::fv::artificialViscosityUSource::addSup
(
    const volScalarField& rho,
    fvBlockMatrix<vector>& eqn,
    const label fieldI

)
{
    dimensionedScalar viscCoeff
    (
        "visc",
        rho.dimensions()*eqn.psi().dimensions()*dimLength*dimLength/dimTime,
        viscosityCoeff_
    );

    volScalarField viscosity
    (
        IOobject
        (
            "artificialViscosity_" + name(),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        viscCoeff,
        zeroGradientFvPatchField<scalar>::typeName
    );
    viscosity.primitiveFieldRef() *= sensor_().valueField();
    viscosity.primitiveFieldRef() =
        min(viscosity.primitiveField(), maxViscosity_);
    viscosity.correctBoundaryConditions();

    if (mesh_.time().outputTime())
    {
        sensor_().writeField();

        if (writeViscosity_ == true)
        {
            viscosity.write();
        }
    }

    eqn +=
        fvm::laplacian(viscosity, eqn.psi())
      + fvc::div
        (
            viscosity*dev2(T(fvc::grad(eqn.psi()))),
            "div((nu*dev2(T(grad(U)))))"
        );
}


void Foam::fv::artificialViscosityUSource::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


// ************************************************************************* //
