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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "turbulenceControlledDecaySource.H"
#include "fvMatrices/fvMatrices.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"
#include "turbulenceModel.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(turbulenceControlledDecaySource, 0);

    addToRunTimeSelectionTable
    (
        option,
        turbulenceControlledDecaySource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::turbulenceControlledDecaySource::turbulenceControlledDecaySource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(sourceName, modelType, dict, obr),
    kAmb_(coeffs_.lookup("kAmb")),
    tAmb_(coeffs_.lookup("tAmb"))
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::turbulenceControlledDecaySource::sourceFields
(
    wordList& fieldNames
)
{
    fieldNames = coeffs_.lookup<wordList>("fields");
}


void Foam::fv::turbulenceControlledDecaySource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (eqn.psi().name() == word(kName_))
    {
        // for k
        const scalar betaStar(0.09);

        eqn += betaStar*tAmb_*kAmb_;
    }

    if (eqn.psi().name() == word(tName_))
    {
        // for omega
        const surfaceScalarField& phi =
            obr_.lookupObject<surfaceScalarField>("phi");

        tmp<volScalarField::Internal> dSource
        (
            new volScalarField::Internal
            (
                IOobject
                (
                    "dSource",
                    mesh().time().timeName(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar
                (
                    "dSource", dimless, 0.0
                )
            )
        );

        tmp<volScalarField> dSourceVsf
        (
            new volScalarField
            (
                IOobject
                (
                    "dSourceVsf",
                    mesh().time().timeName(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar
                (
                    "dSourceVsf", dimless, 0.0
                )
            )
        );

        if (phi.dimensions() == dimVelocity*dimArea)
        {
            const incompressible::turbulenceModel& turbModel =
                obr_.lookupObject<incompressible::turbulenceModel>
                (
                    turbulenceModel::propertiesName
                );

            dSource.ref() = turbModel.decaySource();
        }
        else
        {
            const compressible::turbulenceModel& turbModel =
                obr_.lookupObject<compressible::turbulenceModel>
                (
                    turbulenceModel::propertiesName
                );

            dSource.ref() = turbModel.decaySource();
        }

        dSourceVsf.ref().primitiveFieldRef() = dSource;

        eqn += dSourceVsf*tAmb_*tAmb_;
    }

    if (debug)
    {
        Info<< "calling 2-arg addSup" << endl;
        Info<<"Inserting controlled decay source term for " << eqn.psi().name() << endl;
    }
}


void Foam::fv::turbulenceControlledDecaySource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (eqn.psi().name() == word(kName_))
    {
        // for k
        const scalar betaStar(0.09);
        eqn += rho*betaStar*tAmb_*kAmb_;
    }

    if (eqn.psi().name() == word(tName_))
    {
        // for omega
        const surfaceScalarField& phi =
            obr_.lookupObject<surfaceScalarField>("phi");

        tmp<volScalarField::Internal> dSource
        (
            new volScalarField::Internal
            (
                IOobject
                (
                    "dSource",
                    mesh().time().timeName(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar
                (
                    "dSource", dimless, 0.0
                )
            )
        );

        tmp<volScalarField> dSourceVsf
        (
            new volScalarField
            (
                IOobject
                (
                    "dSourceVsf",
                    mesh().time().timeName(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar
                (
                    "dSourceVsf", dimless, 0.0
                )
            )
        );

        if (phi.dimensions() == dimVelocity*dimArea)
        {
            const incompressible::turbulenceModel& turbModel =
                obr_.lookupObject<incompressible::turbulenceModel>
                (
                    turbulenceModel::propertiesName
                );

            dSource.ref() = turbModel.decaySource();
        }
        else
        {
            const compressible::turbulenceModel& turbModel =
                obr_.lookupObject<compressible::turbulenceModel>
                (
                    turbulenceModel::propertiesName
                );

            dSource.ref() = turbModel.decaySource();
        }

        dSourceVsf.ref().primitiveFieldRef() = dSource;
        eqn += rho*dSourceVsf*tAmb_*tAmb_;
    }

    if (debug)
    {
        Info<< "calling 3-arg addSup" << endl;
        Info<<"Inserting controlled decay source term for " << eqn.psi().name() << endl;
    }
}


bool Foam::fv::turbulenceControlledDecaySource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        kName_ = coeffs_.lookupOrDefault<word>("kName", "k");
        tName_ = coeffs_.lookupOrDefault<word>("tName", "omega");

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
