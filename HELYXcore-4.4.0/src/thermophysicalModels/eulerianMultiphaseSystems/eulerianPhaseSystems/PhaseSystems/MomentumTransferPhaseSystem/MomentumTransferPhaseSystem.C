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

\*---------------------------------------------------------------------------*/

#include "MomentumTransferPhaseSystem.H"

#include "../../EulerianBlendedInterfacialModel/EulerianBlendedInterfacialModel.H"
#include "../../interfacialModels/eulerianDragModels/eulerianDragModel/eulerianDragModel.H"
#include "../../interfacialModels/eulerianVirtualMassModels/eulerianVirtualMassModel/eulerianVirtualMassModel.H"
#include "../../interfacialModels/eulerianLiftModels/eulerianLiftModel/eulerianLiftModel.H"
#include "../../interfacialModels/eulerianWallLubricationModels/eulerianWallLubricationModel/eulerianWallLubricationModel.H"
#include "../../interfacialModels/eulerianTurbulentDispersionModels/eulerianTurbulentDispersionModel/eulerianTurbulentDispersionModel.H"

#include "containers/HashTables/HashPtrTable/HashPtrTable.H"

#include "finiteVolume/fvm/fvmDdt.H"
#include "finiteVolume/fvm/fvmDiv.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "finiteVolume/fvc/fvcDiv.H"
#include "finiteVolume/fvc/fvcSnGrad.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::
MomentumTransferPhaseSystem
(
    const fvMesh& mesh,
    const dictionary& coeffs
)
:
    BasePhaseSystem(mesh, coeffs)
{
    this->createSubModels(dragModels_);
    this->createSubModels(virtualMassModels_);
    this->createSubModels(liftModels_);
    this->createSubModels(wallLubricationModels_);
    this->createSubModels(turbulentDispersionModels_);

    forAllConstIter
    (
        dragModelTable,
        dragModels_,
        dragModelIter
    )
    {
        const eulerianPhasePair& pair(this->phasePairs_[dragModelIter.key()]);

        Kds_.insert
        (
            pair,
            new volScalarField
            (
                IOobject::groupName("Kd", pair.name()),
                dragModelIter()->K()
            )
        );
    }

    forAllConstIter
    (
        virtualMassModelTable,
        virtualMassModels_,
        virtualMassModelIter
    )
    {
        const eulerianPhasePair& pair(this->phasePairs_[virtualMassModelIter.key()]);

        Vms_.insert
        (
            pair,
            new volScalarField
            (
                IOobject::groupName("Vm", pair.name()),
                virtualMassModelIter()->K()
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::
~MomentumTransferPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::Kd
(
    const phasePairKey& key
) const
{
    return dragModels_[key]->K();
}


template<class BasePhaseSystem>
Foam::tmp<Foam::surfaceScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::Kdf
(
    const phasePairKey& key
) const
{
    return dragModels_[key]->Kf();
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::Kd
(
    const Foam::eulerianPhaseModel& phase
) const
{
    tmp<volScalarField> tKd
    (
        volScalarField::New
        (
            IOobject::groupName("Kd", phase.name()),
            this->mesh_,
            dimensionedScalar(dimensionSet(1, -3, -1, 0, 0), 0)
        )
    );

    forAllConstIter
    (
        eulerianPhaseSystem::KdTable,
        Kds_,
        KdIter
    )
    {
        const volScalarField& K(*KdIter());

        const eulerianPhasePair& pair(this->phasePairs_[KdIter.key()]);

        const eulerianPhaseModel* phase1 = &pair.phase1();
        const eulerianPhaseModel* phase2 = &pair.phase2();

        forAllConstIter(eulerianPhasePair, pair, iter)
        {
            if (phase1 == &phase)
            {
                tKd.ref() += K;
            }

            Swap(phase1, phase2);
        }
    }

    return tKd;
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::Vm
(
    const phasePairKey& key
) const
{
    if (virtualMassModels_.found(key))
    {
        return virtualMassModels_[key]->K();
    }
    else
    {
        return volScalarField::New
        (
            eulerianVirtualMassModel::typeName + "-K",
            this->mesh_,
            dimensionedScalar(eulerianVirtualMassModel::dimK, 0)
        );
    }
}


template<class BasePhaseSystem>
Foam::tmp<Foam::surfaceScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::Vmf
(
    const phasePairKey& key
) const
{
    if (virtualMassModels_.found(key))
    {
        return virtualMassModels_[key]->Kf();
    }
    else
    {
        return surfaceScalarField::New
        (
            eulerianVirtualMassModel::typeName + "-Kf",
            this->mesh_,
            dimensionedScalar(eulerianVirtualMassModel::dimK, 0)
        );
    }
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volVectorField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::F
(
    const phasePairKey& key
) const
{
    if (liftModels_.found(key) && wallLubricationModels_.found(key))
    {
        return
            liftModels_[key]->template F<vector>()
          + wallLubricationModels_[key]->template F<vector>();
    }
    else if (liftModels_.found(key))
    {
        return liftModels_[key]->template F<vector>();
    }
    else if (wallLubricationModels_.found(key))
    {
        return wallLubricationModels_[key]->template F<vector>();
    }
    else
    {
        return volVectorField::New
        (
            eulerianLiftModel::typeName + "-F",
            this->mesh_,
            dimensionedVector(eulerianLiftModel::dimF, Zero)
        );
    }
}


template<class BasePhaseSystem>
Foam::tmp<Foam::surfaceScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::Ff
(
    const phasePairKey& key
) const
{
    if (liftModels_.found(key) && wallLubricationModels_.found(key))
    {
        return
            liftModels_[key]->Ff()
          + wallLubricationModels_[key]->Ff();
    }
    else if (liftModels_.found(key))
    {
        return liftModels_[key]->Ff();
    }
    else if (wallLubricationModels_.found(key))
    {
        return wallLubricationModels_[key]->Ff();
    }
    else
    {
        tmp<surfaceScalarField> tFf
        (
            surfaceScalarField::New
            (
                eulerianLiftModel::typeName + "-Ff",
                this->mesh_,
                dimensionedScalar(eulerianLiftModel::dimF*dimArea, 0)
            )
        );

        tFf.ref().setOriented();

        return tFf;
    }
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::D
(
    const phasePairKey& key
) const
{
    if (turbulentDispersionModels_.found(key))
    {
        return turbulentDispersionModels_[key]->D();
    }
    else
    {
        return volScalarField::New
        (
            eulerianTurbulentDispersionModel::typeName + "-D",
            this->mesh_,
            dimensionedScalar(eulerianTurbulentDispersionModel::dimD, 0)
        );
    }
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::eulerianPhaseSystem::momentumTransferTable>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::momentumTransfer() const
{
    // Create a momentum transfer matrix for each phase
    autoPtr<eulerianPhaseSystem::momentumTransferTable> eqnsPtr
    (
        new eulerianPhaseSystem::momentumTransferTable()
    );

    eulerianPhaseSystem::momentumTransferTable& eqns = eqnsPtr();

    forAll(this->phaseModels_, phasei)
    {
        const eulerianPhaseModel& phase = this->phaseModels_[phasei];

        eqns.insert
        (
            phase.name(),
            new fvVectorMatrix(phase.U(), dimMass*dimVelocity/dimTime)
        );
    }

    // Update the drag coefficients
    forAllConstIter
    (
        dragModelTable,
        dragModels_,
        dragModelIter
    )
    {
        *Kds_[dragModelIter.key()] = dragModelIter()->K();
    }

    // Add the implicit part of the drag force
    forAllConstIter
    (
        eulerianPhaseSystem::KdTable,
        Kds_,
        KdIter
    )
    {
        const volScalarField& K(*KdIter());

        const eulerianPhasePair& pair(this->phasePairs_[KdIter.key()]);

        const eulerianPhaseModel* phase = &pair.phase1();
        const eulerianPhaseModel* otherPhase = &pair.phase2();

        forAllConstIter(eulerianPhasePair, pair, iter)
        {
            tmp<volVectorField> tU(phase->U());
            const volVectorField& U = tU();

            *eqns[phase->name()] -= fvm::Sp(K, U);

            Swap(phase, otherPhase);
        }
    }

    // Update the virtual mass coefficients
    forAllConstIter
    (
        virtualMassModelTable,
        virtualMassModels_,
        virtualMassModelIter
    )
    {
        *Vms_[virtualMassModelIter.key()] = virtualMassModelIter()->K();
    }

    // Add the virtual mass force
    forAllConstIter
    (
        eulerianPhaseSystem::VmTable,
        Vms_,
        VmIter
    )
    {
        const volScalarField& Vm(*VmIter());

        const eulerianPhasePair& pair(this->phasePairs_[VmIter.key()]);

        const eulerianPhaseModel* phase = &pair.phase1();
        const eulerianPhaseModel* otherPhase = &pair.phase2();

        forAllConstIter(eulerianPhasePair, pair, iter)
        {
            tmp<volVectorField> tU(phase->U());
            const volVectorField& U = tU();

            tmp<surfaceScalarField> tphiv(phase->phiv());
            const surfaceScalarField& phiv = tphiv();

            *eqns[phase->name()] -=
                Vm
               *(
                    fvm::ddt(U)
                  + fvm::div(phiv, U)
                  - fvm::Sp(fvc::div(phiv), U)
                  - otherPhase->DUDt()
                )
              + this->fvOptions().MRFDDt(Vm, U - otherPhase->U());

            Swap(phase, otherPhase);
        }
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::volVectorField& Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::setF
(
    PtrList<volVectorField>& Fs, const label phasei
) const
{
    if (!Fs.set(phasei))
    {
        Fs.set
        (
            phasei,
            new volVectorField
            (
                IOobject
                (
                    eulerianLiftModel::typeName + "-F",
                    this->mesh_.time().timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                this->mesh_,
                dimensionedVector(eulerianLiftModel::dimF, Zero)
            )
        );
    }

    return Fs[phasei];
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::PtrList<Foam::volVectorField>>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::Fs() const
{
    autoPtr<PtrList<volVectorField>> tFs
    (
        new PtrList<volVectorField>(this->phases().size())
    );
    PtrList<volVectorField>& Fs = tFs();

    // Add the lift force
    forAllConstIter
    (
        liftModelTable,
        liftModels_,
        liftModelIter
    )
    {
        const volVectorField F(liftModelIter()->F<vector>());

        const eulerianPhasePair& pair(this->phasePairs_[liftModelIter.key()]);

        setF(Fs, pair.phase1().index()) += F;
        setF(Fs, pair.phase2().index()) -= F;
    }

    // Add the wall lubrication force
    forAllConstIter
    (
        wallLubricationModelTable,
        wallLubricationModels_,
        wallLubricationModelIter
    )
    {
        const volVectorField F(wallLubricationModelIter()->F<vector>());

        const eulerianPhasePair&
            pair(this->phasePairs_[wallLubricationModelIter.key()]);

        setF(Fs, pair.phase1().index()) += F;
        setF(Fs, pair.phase2().index()) -= F;
    }

    return tFs;
}


template<class BasePhaseSystem>
Foam::surfaceScalarField&
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::setPhiD
(
    PtrList<surfaceScalarField>& phiDs, const label phasei
) const
{
    if (!phiDs.set(phasei))
    {
        phiDs.set
        (
            phasei,
            new surfaceScalarField
            (
                IOobject
                (
                    eulerianTurbulentDispersionModel::typeName + "-phiD",
                    this->mesh_.time().timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                this->mesh_,
                dimensionedScalar
                (
                    dimTime*dimArea*eulerianTurbulentDispersionModel::dimF/dimDensity,
                    0
                )
            )
        );

        phiDs[phasei].setOriented();
    }

    return phiDs[phasei];
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::PtrList<Foam::surfaceScalarField>>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::phiDs
(
    const PtrList<volScalarField>& rAUs
) const
{
    autoPtr<PtrList<surfaceScalarField>> tphiDs
    (
        new PtrList<surfaceScalarField>(this->phases().size())
    );
    PtrList<surfaceScalarField>& phiDs = tphiDs();

    // Add the turbulent dispersion force
    forAllConstIter
    (
        turbulentDispersionModelTable,
        turbulentDispersionModels_,
        turbulentDispersionModelIter
    )
    {
        const eulerianPhasePair&
            pair(this->phasePairs_[turbulentDispersionModelIter.key()]);

        const volScalarField D(turbulentDispersionModelIter()->D());
        const surfaceScalarField snGradAlpha1
        (
            fvc::snGrad(pair.phase1().volFrac())*this->mesh_.magSf()
        );

        setPhiD(phiDs, pair.phase1().index()) +=
            fvc::interpolate(rAUs[pair.phase1().index()]*D)*snGradAlpha1;
        setPhiD(phiDs, pair.phase2().index()) -=
            fvc::interpolate(rAUs[pair.phase2().index()]*D)*snGradAlpha1;
    }

    return tphiDs;
}


template<class BasePhaseSystem>
bool Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::read()
{
    if (BasePhaseSystem::read())
    {
        bool readOK = true;

        // Read models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
