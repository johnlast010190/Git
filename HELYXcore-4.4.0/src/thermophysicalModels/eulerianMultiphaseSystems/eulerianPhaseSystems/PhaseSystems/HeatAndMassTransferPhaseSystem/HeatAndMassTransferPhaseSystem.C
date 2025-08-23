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

#include "HeatAndMassTransferPhaseSystem.H"

#include "../../EulerianBlendedInterfacialModel/EulerianBlendedInterfacialModel.H"
#include "../../interfacialModels/eulerianHeatTransferModels/eulerianHeatTransferModel/eulerianHeatTransferModel.H"
#include "../../interfacialCompositionModels/eulerianMassTransferModels/eulerianMassTransferModel/eulerianMassTransferModel.H"

#include "containers/HashTables/HashPtrTable/HashPtrTable.H"

#include "finiteVolume/fvc/fvcDiv.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::
HeatAndMassTransferPhaseSystem
(
    const fvMesh& mesh,
    const dictionary& coeffs
)
:
    BasePhaseSystem(mesh, coeffs)
{
    this->createSubModels(heatTransferModels_);
    this->createSubModels(massTransferModels_);

    forAllConstIter
    (
        eulerianPhaseSystem::phasePairTable,
        this->phasePairs_,
        phasePairIter
    )
    {
        const eulerianPhasePair& pair(phasePairIter());

        if (pair.ordered())
        {
            continue;
        }

        // Initialy assume no mass transfer

        dmdt_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("dmdt", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimDensity/dimTime, 0)
            )
        );

        dmdtExplicit_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("dmdtExplicit", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh()
                ),
                this->mesh(),
                dimensionedScalar(dimDensity/dimTime, 0)
            )
        );

        volScalarField H1(heatTransferModels_[pair][pair.first()]->K());
        volScalarField H2(heatTransferModels_[pair][pair.second()]->K());

        Tf_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("Tf", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                (
                    H1*pair.phase1().thermo().T()
                  + H2*pair.phase2().thermo().T()
                )
               /max
                (
                    H1 + H2,
                    dimensionedScalar(eulerianHeatTransferModel::dimK, SMALL)
                ),
                zeroGradientFvPatchScalarField::typeName
            )
        );
        Tf_[pair]->correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::
~HeatAndMassTransferPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
bool Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::transfersMass
(
    const eulerianPhaseModel& phase
) const
{
    return true;
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::dmdt
(
    const phasePairKey& key
) const
{
    const scalar dmdtSign(Pair<word>::compare(dmdt_.find(key).key(), key));

    return dmdtSign**dmdt_[key];
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::dmdt
(
    const Foam::eulerianPhaseModel& phase
) const
{
    tmp<volScalarField> tdmdt
    (
        volScalarField::New
        (
            IOobject::groupName("dmdt", phase.name()),
            this->mesh_,
            dimensionedScalar(dimDensity/dimTime, 0)
        )
    );

    forAllConstIter
    (
        eulerianPhaseSystem::phasePairTable,
        this->phasePairs_,
        phasePairIter
    )
    {
        const eulerianPhasePair& pair(phasePairIter());

        if (pair.ordered())
        {
            continue;
        }

        const eulerianPhaseModel* phase1 = &pair.phase1();
        const eulerianPhaseModel* phase2 = &pair.phase2();

        forAllConstIter(eulerianPhasePair, pair, iter)
        {
            if (phase1 == &phase)
            {
                tdmdt.ref() += this->dmdt(pair);
            }

            Swap(phase1, phase2);
        }
    }

    return tdmdt;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::eulerianPhaseSystem::momentumTransferTable>
Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::momentumTransfer() const
{
    autoPtr<eulerianPhaseSystem::momentumTransferTable>
        eqnsPtr(BasePhaseSystem::momentumTransfer());

    eulerianPhaseSystem::momentumTransferTable& eqns = eqnsPtr();

    // Source term due to mass trasfer
    forAllConstIter
    (
        eulerianPhaseSystem::phasePairTable,
        this->phasePairs_,
        phasePairIter
    )
    {
        const eulerianPhasePair& pair(phasePairIter());

        if (pair.ordered())
        {
            continue;
        }

        tmp<volVectorField> tU1(pair.phase1().U());
        tmp<volVectorField> tU2(pair.phase2().U());
        const volVectorField& U1 = tU1();
        const volVectorField& U2 = tU2();

        const volScalarField dmdt(this->dmdt(pair));
        const volScalarField dmdt21(posPart(dmdt));
        const volScalarField dmdt12(negPart(dmdt));

        *eqns[pair.phase1().name()] += dmdt21*U2 - fvm::Sp(dmdt21, U1);
        *eqns[pair.phase2().name()] -= dmdt12*U1 - fvm::Sp(dmdt12, U2);
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::eulerianPhaseSystem::heatTransferTable>
Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::heatTransfer() const
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

    // Heat transfer with the interface
    forAllConstIter
    (
        heatTransferModelTable,
        heatTransferModels_,
        heatTransferModelIter
    )
    {
        const eulerianPhasePair& pair
        (
            this->phasePairs_[heatTransferModelIter.key()]
        );

        const eulerianPhaseModel* phase = &pair.phase1();
        const eulerianPhaseModel* otherPhase = &pair.phase2();

        const volScalarField& Tf(*Tf_[pair]);

        const volScalarField K1
        (
            heatTransferModelIter()[pair.first()]->K()
        );
        const volScalarField K2
        (
            heatTransferModelIter()[pair.second()]->K()
        );
        const volScalarField KEff
        (
            K1*K2
           /max
            (
                K1 + K2,
                dimensionedScalar(eulerianHeatTransferModel::dimK, SMALL)
            )
        );

        const volScalarField* K = &K1;
        const volScalarField* otherK = &K2;

        forAllConstIter(eulerianPhasePair, pair, iter)
        {
            const volScalarField& he = phase->thermo().he();
            const volScalarField& Cpv = phase->thermo().Cpv();

            *eqns[phase->name()] +=
                (*K)*(Tf - phase->thermo().T())
              + KEff/Cpv*he - fvm::Sp(KEff/Cpv, he);

            Swap(phase, otherPhase);
            Swap(K, otherK);
        }
    }

    // Source term due to mass transfer
    forAllConstIter
    (
        eulerianPhaseSystem::phasePairTable,
        this->phasePairs_,
        phasePairIter
    )
    {
        const eulerianPhasePair& pair(phasePairIter());

        if (pair.ordered())
        {
            continue;
        }

        const eulerianPhaseModel& phase1 = pair.phase1();
        const eulerianPhaseModel& phase2 = pair.phase2();

        const volScalarField& he1(phase1.thermo().he());
        const volScalarField& he2(phase2.thermo().he());

        const volScalarField& K1(phase1.K());
        const volScalarField& K2(phase2.K());

        const volScalarField dmdt(this->dmdt(pair));
        const volScalarField dmdt21(posPart(dmdt));
        const volScalarField dmdt12(negPart(dmdt));
        const volScalarField& Tf(*Tf_[pair]);

        *eqns[phase1.name()] +=
            dmdt21*(phase1.thermo().he(phase1.thermo().p(), Tf))
          - fvm::Sp(dmdt21, he1)
          + dmdt21*(K2 - K1);

        *eqns[phase2.name()] -=
            dmdt12*(phase2.thermo().he(phase2.thermo().p(), Tf))
          - fvm::Sp(dmdt12, he2)
          + dmdt12*(K1 - K2);
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
bool Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::read()
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
