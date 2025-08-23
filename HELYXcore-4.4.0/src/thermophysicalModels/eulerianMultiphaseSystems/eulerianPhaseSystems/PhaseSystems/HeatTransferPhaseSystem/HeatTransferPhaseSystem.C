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
    (c) 2015-2022 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "HeatTransferPhaseSystem.H"

#include "../../EulerianBlendedInterfacialModel/EulerianBlendedInterfacialModel.H"

// Must explicitly specify that we want eulerianHeatTransferModel from eulerianPhaseSystem, not
// multiphase system!
#include "../../interfacialModels/eulerianHeatTransferModels/eulerianHeatTransferModel/eulerianHeatTransferModel.H"

#include "containers/HashTables/HashPtrTable/HashPtrTable.H"

#include "finiteVolume/fvc/fvcDiv.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::HeatTransferPhaseSystem<BasePhaseSystem>::HeatTransferPhaseSystem
(
    const fvMesh& mesh,
    const dictionary& coeffs
)
:
    BasePhaseSystem(mesh, coeffs)
{
    this->createSubModels(heatTransferModels_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::HeatTransferPhaseSystem<BasePhaseSystem>::~HeatTransferPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
bool Foam::HeatTransferPhaseSystem<BasePhaseSystem>::transfersMass
(
    const eulerianPhaseModel& phase
) const
{
    return false;
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::HeatTransferPhaseSystem<BasePhaseSystem>::dmdt
(
    const phasePairKey& key
) const
{
    return volScalarField::New
    (
        IOobject::groupName("dmdt", this->phasePairs_[key]->name()),
        this->mesh(),
        dimensionedScalar(dimDensity/dimTime, 0)
    );
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::HeatTransferPhaseSystem<BasePhaseSystem>::dmdt
(
    const Foam::eulerianPhaseModel& phase
) const
{
    return tmp<volScalarField>(nullptr);
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::eulerianPhaseSystem::heatTransferTable>
Foam::HeatTransferPhaseSystem<BasePhaseSystem>::heatTransfer() const
{
    autoPtr<eulerianPhaseSystem::heatTransferTable> eqnsPtr
    (
        new eulerianPhaseSystem::heatTransferTable()
    );

    eulerianPhaseSystem::heatTransferTable& eqns = eqnsPtr();

    forAll(this->phaseModels_, phasei)
    {
        const eulerianPhaseModel& phase = this->phaseModels_[phasei];

        eqns.insert
        (
            phase.name(),
            new fvScalarMatrix(phase.thermo().he(), dimEnergy/dimTime)
        );
    }

    forAllConstIter
    (
        heatTransferModelTable,
        heatTransferModels_,
        heatTransferModelIter
    )
    {
        const volScalarField K(heatTransferModelIter()->K());

        const eulerianPhasePair& pair(this->phasePairs_[heatTransferModelIter.key()]);

        const eulerianPhaseModel* phase = &pair.phase1();
        const eulerianPhaseModel* otherPhase = &pair.phase2();

        forAllConstIter(eulerianPhasePair, pair, iter)
        {
            const volScalarField& he = phase->thermo().he();
            const volScalarField& Cpv = phase->thermo().Cpv();

            *eqns[phase->name()] +=
                K*(otherPhase->thermo().T() - phase->thermo().T() + he/Cpv)
              - fvm::Sp(K/Cpv, he);

            Swap(phase, otherPhase);
        }
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::eulerianPhaseSystem::massTransferTable>
Foam::HeatTransferPhaseSystem<BasePhaseSystem>::massTransfer() const
{
    autoPtr<eulerianPhaseSystem::massTransferTable> eqnsPtr
    (
        new eulerianPhaseSystem::massTransferTable()
    );

    return eqnsPtr;
}


template<class BasePhaseSystem>
bool Foam::HeatTransferPhaseSystem<BasePhaseSystem>::read()
{
    if (BasePhaseSystem::read())
    {
        bool readOK = true;

        // Models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
