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
    (c) 2023-2024 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "eulerianPhaseSystem.H"
#include "../interfacialModels/eulerianAspectRatioModels/eulerianAspectRatioModel/eulerianAspectRatioModel.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "finiteVolume/fvc/fvcDdt.H"
#include "eulerianPhaseSystems/eulerianPhaseModel/MovingPhaseModel/phaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(eulerianPhaseSystem, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::eulerianPhaseSystem::calcPhiv
(
    const phaseModelList& eulerianPhaseModels
) const
{
    tmp<surfaceScalarField> tmpPhiv
    (
        new surfaceScalarField
        (
            "phiv",
            fvc::interpolate(eulerianPhaseModels[0].volFrac())*eulerianPhaseModels[0].phiv()
        )
    );

    for (label phasei=1; phasei<eulerianPhaseModels.size(); phasei++)
    {
        tmpPhiv.ref() +=
            fvc::interpolate(eulerianPhaseModels[phasei].volFrac())*eulerianPhaseModels[phasei].phiv();
    }

    return tmpPhiv;
}


void Foam::eulerianPhaseSystem::generatePairs()
{
    wordList phases(thermo_.phases());
    forAll(phases, phasei)
    {
        forAll(phases, phasej)
        {
            if (phasei == phasej)
            {
                continue;
            }
            phasePairKey key
            (
                phases[phasei], phases[phasej], false
            );
            phasePairKey orderedKey
            (
                phases[phasei], phases[phasej], true
            );
            if (phasei < phasej)
            {
                phasePairs_.insert
                (
                    key,
                    autoPtr<eulerianPhasePair>
                    (
                        new eulerianPhasePair
                        (
                            phaseModels_[key.first()],
                            phaseModels_[key.second()]
                        )
                    )
                );
            }
            phasePairs_.insert
            (
                orderedKey,
                autoPtr<eulerianPhasePair>
                (
                    new orderedEulerianPhasePair
                    (
                        phaseModels_[orderedKey.first()],
                        phaseModels_[orderedKey.second()]
                    )
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianPhaseSystem::eulerianPhaseSystem
(
    const fvMesh& mesh,
    const dictionary& coeffs
)
:
    regIOobject
    (
        IOobject
        (
            typeName,
            mesh.time().constant(),
            mesh
        )
    ),
    mesh_(mesh),
    thermo_(refCast<multiphaseThermo>(multiphaseThermo::lookupOrCreate(mesh_))),
    phaseModels_(thermo_.lookup("phases"), eulerianPhaseModel::iNew(*this)),
    phiv_(calcPhiv(phaseModels_)),
    dpdt_
    (
        IOobject
        (
            "dpdt",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimPressure/dimTime, 0)
    )
{
    phiv_.writeOpt() = IOobject::AUTO_WRITE;

    // Blending methods
    forAllConstIter(dictionary, coeffs.subDict("blending"), iter)
    {
        blendingMethods_.insert
        (
            iter().dict().dictName(),
            eulerianBlendingMethod::New
            (
                iter().dict(),
                phaseModels_.toc()
            )
        );
    }

    generatePairs();

    // Sub-models
    createSubModels(aspectRatioModels_, true);

    correctKinematics();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianPhaseSystem::~eulerianPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::eulerianPhaseSystem::rho() const
{
    tmp<volScalarField> tmpRho
    (
        phaseModels_[0].volFrac()*phaseModels_[0].rho()
    );

    for (label phasei=1; phasei<phaseModels_.size(); phasei++)
    {
        tmpRho.ref() +=
            phaseModels_[phasei].volFrac()*phaseModels_[phasei].rho();
    }

    return tmpRho;
}


Foam::tmp<Foam::volVectorField> Foam::eulerianPhaseSystem::U() const
{
    tmp<volVectorField> tmpU
    (
        phaseModels_[0].volFrac()*phaseModels_[0].U()
    );

    for (label phasei=1; phasei<phaseModels_.size(); phasei++)
    {
        tmpU.ref() += phaseModels_[phasei].volFrac()*phaseModels_[phasei].U();
    }

    return tmpU;
}


Foam::tmp<Foam::volScalarField> Foam::eulerianPhaseSystem::muEff() const
{
    tmp<volScalarField> tmpmuEff
    (
        phaseModels_[0].volFrac()*phaseModels_[0].turbulence().muEff()
    );

    for (label phasei=1; phasei<phaseModels_.size(); phasei++)
    {
        tmpmuEff.ref() +=
            phaseModels_[phasei].volFrac()
           *phaseModels_[phasei].turbulence().muEff();
    }

    return tmpmuEff;
}


Foam::tmp<Foam::volScalarField>
Foam::eulerianPhaseSystem::E(const phasePairKey& key) const
{
    if (aspectRatioModels_.found(key))
    {
        return aspectRatioModels_[key]->E();
    }
    else
    {
        return volScalarField::New
        (
            eulerianAspectRatioModel::typeName + "-E",
            this->mesh_,
            dimensionedScalar(dimless, 1)
        );
    }
}


void Foam::eulerianPhaseSystem::solve()
{}


void Foam::eulerianPhaseSystem::correct()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correct();
    }
}


void Foam::eulerianPhaseSystem::correctKinematics()
{
    bool updateDpdt = false;

    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctKinematics();

        updateDpdt = updateDpdt || phaseModels_[phasei].thermo().dpdt();
    }

    // Update the pressure time-derivative if required
    if (updateDpdt)
    {
        dpdt_ = fvc::ddt(phaseModels_.begin()().thermo().p());
    }
}


void Foam::eulerianPhaseSystem::correctThermo()
{
    thermo_.correct();
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctThermo();
    }
}


void Foam::eulerianPhaseSystem::correctTurbulence()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctTurbulence();
    }
}


void Foam::eulerianPhaseSystem::correctEnergyTransport()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctEnergyTransport();
    }
}


bool Foam::eulerianPhaseSystem::read()
{
    bool readOK = true;

    forAll(phaseModels_, phasei)
    {
        readOK &= phaseModels_[phasei].read();
    }

    // models ...

    return readOK;
}


// ************************************************************************* //
