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
    (c) 2022-2023 OpenFOAM Foundation
    (c) 2024 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "phaseTurbulenceStabilisation.H"
#include "eulerianPhaseSystems/eulerianPhaseSystem/eulerianPhaseSystem.H"
#include "turbulenceModel.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "eulerianPhaseSystems/eulerianPhaseModel/MovingPhaseModel/phaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(phaseTurbulenceStabilisation, 0);

        addToRunTimeSelectionTable
        (
            option,
            phaseTurbulenceStabilisation,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::phaseTurbulenceStabilisation::addAlphaRhoSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volScalarField& field,
    fvMatrix<scalar>& eqn,
    tmp<volScalarField> (turbulenceModel::*psi)() const
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    const fvMesh& mesh = this->mesh();

    //TODO: In future, if stationary phases are supported (porosity?) this
    // should only be the moving phases. At the moment it is all phases
    const eulerianPhaseSystem::phaseModelList& movingPhases =
        phaseSystem_->phases();

    volScalarField::Internal transferRate
    (
        tmp<volScalarField::Internal>
        (
            new volScalarField::Internal
            (
                IOobject
                (
                    "transferRate",
                    mesh_.time().timeName(),
                    eqn.psi().db()
                ),
                mesh,
                dimensionedScalar(dimless/dimTime, 0)
            )
        )
    );
    volScalarField::Internal psiTransferRate
    (
        tmp<volScalarField::Internal>
        (
            new volScalarField::Internal
            (
                IOobject
                (
                    "psiTransferRate",
                    mesh_.time().timeName(),
                    eqn.psi().db()
                ),
                mesh,
                dimensionedScalar(eqn.psi().dimensions()/dimTime, 0)
            )
        )
    );

    bool otherTurbModelFound = false;
    forAll(movingPhases, phasei)
    {
        if (movingPhases[phasei].name() != phaseName_)
        {
            const turbulenceModel* turbulence =
                obr_.lookupObjectPtr<turbulenceModel>
                (
                    IOobject::groupName
                    (
                        turbulenceModel::propertiesName,
                        movingPhases[phasei].name()
                    )
                );

            if
            (
                turbulence &&
               !isA<laminareulerianPhaseModelPhaseCompressibleTurbulenceModel>
                (
                    *turbulence
                )
            )
            {
                const volScalarField::Internal phaseTransferRate
                (
                    min
                    (
                        turbulence->epsilon()/turbulence->k(),
                        1.0/mesh_.time().deltaT()
                    )*movingPhases[phasei].volFrac()
                );

                transferRate += phaseTransferRate;
                psiTransferRate += phaseTransferRate*(turbulence->*psi)()();
                otherTurbModelFound = true;
            }
        }
    }
    if (!otherTurbModelFound)
    {
        WarningInFunction
            << "No turbulence model found in other phases to use for "
            << "turbulence stabilisation of phase '" << phaseName_ << "'"
            << nl << endl;
    }

    const volScalarField::Internal transferCoeff
    (
        max(alphaInversion_ - alpha(), scalar(0))*rho()
    );

    eqn += transferCoeff*psiTransferRate;
    eqn -= fvm::Sp(transferCoeff*transferRate, eqn.psi());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::phaseTurbulenceStabilisation::phaseTurbulenceStabilisation
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    option(sourceName, modelType, dict, obr),
    phaseName_(dict.lookup("phase")),
    alphaInversion_("alphaInversion", dimless, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::phaseTurbulenceStabilisation::sourceFields(wordList& fieldNames)
{
    phaseSystem_ =
        &obr_.lookupObject<eulerianPhaseSystem>(eulerianPhaseSystem::typeName);
    turbulence_ =
        &obr_.lookupObject<turbulenceModel>
        (
            IOobject::groupName(turbulenceModel::propertiesName, phaseName_)
        );

    const word kName(IOobject::groupName("k", phaseName_));
    const word epsilonName(IOobject::groupName("epsilon", phaseName_));
    const word omegaName(IOobject::groupName("omega", phaseName_));

    if (obr_.foundObject<volScalarField>(kName))
    {
        fieldNames.append(kName);
    }

    if (obr_.foundObject<volScalarField>(epsilonName))
    {
        fieldNames.append(epsilonName);
    }

    if (obr_.foundObject<volScalarField>(omegaName))
    {
        fieldNames.append(omegaName);
    }
}


void Foam::fv::phaseTurbulenceStabilisation::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    const volScalarField& field(eqn.psi());
    if (field.name() == IOobject::groupName("k", phaseName_))
    {
        addAlphaRhoSup
        (
            alpha,
            rho,
            field,
            eqn,
            &turbulenceModel::k
        );
    }
    else if (field.name() == IOobject::groupName("epsilon", phaseName_))
    {
        addAlphaRhoSup
        (
            alpha,
            rho,
            field,
            eqn,
            &turbulenceModel::epsilon
        );
    }
    else if (field.name() == IOobject::groupName("omega", phaseName_))
    {
        addAlphaRhoSup
        (
            alpha,
            rho,
            field,
            eqn,
            &turbulenceModel::omega
        );
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << field.name() << " is not implemented"
            << exit(FatalError);
    }
}


// ************************************************************************* //
