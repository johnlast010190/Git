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
    (c) 2015-2017 OpenFOAM Foundation
    (c) 2024 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "../EulerianBlendedInterfacialModel/EulerianBlendedInterfacialModel.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class modelType>
void Foam::eulerianPhaseSystem::createSubModels
(
    HashTable
    <
        autoPtr<modelType>,
        phasePairKey,
        phasePairKey::hash
    >& models,
    const bool ordered,
    const word& phaseName
)
{
    wordList phases(thermo_.phases());
    forAll(phases, phasei)
    {
        forAll(phases, phasej)
        {
            if (phasei == phasej)
            {
                if (ordered)
                {
                    continue;
                }
                else
                {
                    break;
                }
            }

            phasePairKey key
            (
                phases[phasei], phases[phasej], ordered
            );

            word modelTypeName =
                IOobject::groupName
                (
                    modelType::typeName+(ordered ? "" : "Symmetric"),
                    phaseName
                );
            const dictionary* dict =
                thermo_.phasePairDictContaining
                (
                    modelTypeName,
                    phases[phasei],
                    phases[phasej],
                    ordered
                );

            if (dict)
            {
                models.insert
                (
                    key,
                    modelType::New
                    (
                        modelTypeName,
                        dict,
                        phasePairs_[key]
                    )
                );
            }
        }
    }
}


template<class modelType>
void Foam::eulerianPhaseSystem::createSubModels
(
    HashTable
    <
        autoPtr<EulerianBlendedInterfacialModel<modelType>>,
        phasePairKey,
        phasePairKey::hash
    >& models,
    const word& phaseName
)
{
    typedef
        HashTable<autoPtr<modelType>, phasePairKey, phasePairKey::hash>
        modelTypeTable;

    modelTypeTable tempModels;
    createSubModels(tempModels, true, phaseName);
    modelTypeTable tempSymmModels;
    createSubModels(tempSymmModels, false, phaseName);

    const word modelName(IOobject::groupName(modelType::typeName, phaseName));
    const eulerianBlendingMethod& blending
    (
        blendingMethods_.found(modelName)
      ? blendingMethods_[modelName]
      : blendingMethods_["default"]
    );

    autoPtr<modelType> noModel(nullptr);

    forAllConstIter(typename modelTypeTable, tempModels, iter)
    {
        if (!iter().valid())
        {
            continue;
        }

        const phasePairKey key(iter.key().first(), iter.key().second(), false);
        const phasePairKey key1In2(key.first(), key.second(), true);
        const phasePairKey key2In1(key.second(), key.first(), true);

        models.insert
        (
            key,
            autoPtr<EulerianBlendedInterfacialModel<modelType>>
            (
                new EulerianBlendedInterfacialModel<modelType>
                (
                    phaseModels_[key.first()],
                    phaseModels_[key.second()],
                    blending,
                    tempSymmModels.found(key) ? tempSymmModels[key] : noModel,
                    tempModels.found(key1In2) ? tempModels[key1In2] : noModel,
                    tempModels.found(key2In1) ? tempModels[key2In1] : noModel
                )
            )
        );
    }
}


template<class modelType>
void Foam::eulerianPhaseSystem::createSubModels
(
    HashTable
    <
        HashTable<autoPtr<EulerianBlendedInterfacialModel<modelType>>>,
        phasePairKey,
        phasePairKey::hash
    >& models
)
{
    typedef
        HashTable
        <
            autoPtr<EulerianBlendedInterfacialModel<modelType>>,
            phasePairKey,
            phasePairKey::hash
        >
        modelTypeTable;

    forAll(phaseModels_, phasei)
    {
        modelTypeTable tempModels;
        createSubModels
        (
            tempModels,
            phaseModels_[phasei].name()
        );

        forAllConstIter(typename modelTypeTable, tempModels, tempModelIter)
        {
            const phasePairKey key(tempModelIter.key());

            if (!models.found(key))
            {
                models.insert
                (
                    key,
                    HashTable<autoPtr<EulerianBlendedInterfacialModel<modelType>>>()
                );
            }

            models[tempModelIter.key()].insert
            (
                phaseModels_[phasei].name(),
                *tempModelIter
            );
        }
    }
}


template<class modelType>
const modelType& Foam::eulerianPhaseSystem::lookupSubModel(const eulerianPhasePair& key) const
{
    return
        mesh().lookupObject<modelType>
        (
            IOobject::groupName(modelType::typeName, key.name())
        );
}


template<class modelType>
const modelType& Foam::eulerianPhaseSystem::lookupSubModel
(
    const eulerianPhaseModel& dispersed,
    const eulerianPhaseModel& continuous
) const
{
    return lookupSubModel<modelType>(orderedEulerianPhasePair(dispersed, continuous));
}


// ************************************************************************* //
