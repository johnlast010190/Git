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
    (c) 2022-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "itemStorage.H"
#include "itemDataSetProviders/itemDataSetProvider.H"
#include "infos/items/group/boundariesGroupInfo.H"
#include "infos/items/invalid/missingItemInfo.H"
#include "postDict/postDictKeys.H"


namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void ItemStorage::addSourceToProvider(ItemDataSetProvider* provider, const Id &sourceId, const std::string& sceneName)
{
    const ItemInfo* sourceItemInfo = getItemInfo(sourceId, sceneName);
    ItemDataSetProvider* sourceProvider = getOrCreateItemProvider(sourceItemInfo, sceneName);
    provider->addSource(sourceProvider);
}

void ItemStorage::addSourcesToProvider(ItemDataSetProvider* provider, const ItemInfo *itemInfo, const std::string& sceneName)
{
    if (provider->acceptsMeshBoundarySources())
    {
        if (itemInfo->isMeshBoundary())
        {
            FatalErrorInFunction << "Item of type " << typeid(*provider).name() << " both accepts mesh boundaries and is a mesh boundary, which is not permitted" << abort(FatalError);
        }
        const std::vector<Id>& allBoundaries = meshBoundaries_[itemInfo->getId().region].allBoundaries;
        for (const Id &boundaryId : allBoundaries) {
            addSourceToProvider(provider, boundaryId, sceneName);
        }
    }
    else if (provider->acceptsProcessBoundarySources())
    {
        if (itemInfo->isProcessBoundary())
        {
            FatalErrorInFunction << "Item of type " << typeid(*provider).name() << " both accepts process boundaries and is a process boundary, which is not permitted" << abort(FatalError);
        }
        for (const Id &boundaryId : meshBoundaries_[itemInfo->getId().region].processBoundaries) {
            addSourceToProvider(provider, boundaryId, sceneName);
        }
    }
    for (const Id &sourceId : itemInfo->getSources()) {
        addSourceToProvider(provider, sourceId, sceneName);
    }
}

std::vector<ItemDataSetProvider*> ItemStorage::getAllMeshBoundaryProviders(const std::string& region)
{
    std::vector<ItemDataSetProvider*> providers;
    const std::vector<Id>& allBoundaries = meshBoundaries_[region].allBoundaries;
    for (const Id &boundaryId : allBoundaries) {
        const ItemInfo* sourceItemInfo = getBaseItemInfo(boundaryId);
        providers.push_back(getOrCreateItemProvider(sourceItemInfo, ""));
    }
    return providers;
}

std::vector<ItemDataSetProvider*> ItemStorage::getProcessBoundaryProviders(const std::string& region)
{
    std::vector<ItemDataSetProvider*> providers;
    const std::vector<Id>& allBoundaries = meshBoundaries_[region].processBoundaries;
    for (const Id &boundaryId : allBoundaries) {
        const ItemInfo* sourceItemInfo = getBaseItemInfo(boundaryId);
        providers.push_back(getOrCreateItemProvider(sourceItemInfo, ""));
    }
    return providers;
}

ItemDataSetProvider *ItemStorage::getOrCreateItemProvider(const ItemInfo *itemInfo, const std::string &sceneName)
{
    ItemInfoIndex itemInfoIndex(itemInfo, this, sceneName);
    ItemDataSetProvider* provider = getItemProviderIfExists(itemInfoIndex);
    if (!provider)
    {
        provider = itemInfoIndex.itemInfo->createDataSetProvider();
        addSourcesToProvider(provider, itemInfo, sceneName);
        //DebugInFunction << "!!! Created provider for " << itemInfo->getId().name << " - " << sceneName << endl;
        providersFromItemInfo_.emplace(itemInfoIndex, std::unique_ptr<ItemDataSetProvider>(provider));
    }
    return provider;
}

ItemDataSetProvider *ItemStorage::getItemProviderIfExists(const ItemInfoIndex& itemInfoIndex) const
{
    auto iter = providersFromItemInfo_.find(itemInfoIndex);
    if (iter != providersFromItemInfo_.end())
    {
        return iter->second.get();
    }
    else
    {
        return nullptr;
    }
}

ItemDataSetProvider* ItemStorage::createItemProvider(
    const ItemInfo *itemInfo,
    const std::string &sceneName,
    const ItemRequirements& extraRequirements,
    InstanceServices* instanceServices
)
{
    ItemDataSetProvider* provider = getOrCreateItemProvider(itemInfo, sceneName);
    provider->setInstanceServicesForThisAndSources(instanceServices);
    provider->addItemRequirementsForThisAndSources(itemInfo->getItemRequirements());
    provider->addItemRequirementsForThisAndSources(extraRequirements);
    return provider;
}

ItemDataSetProvider* ItemStorage::getItemProvider(const ItemInfo *itemInfo, const std::string& sceneName) const {
    ItemInfoIndex itemInfoIndex(itemInfo, this, sceneName);
    ItemDataSetProvider* provider = getItemProviderIfExists(itemInfoIndex);
    if (!provider)
    {
        FatalErrorInFunction << "Provider for item " << itemInfo->getId().name
                             << " was requested but was not created."
                             << abort(FatalError);
    }
    return provider;
}

ItemDataSetProvider *ItemStorage::getItemProvider(const Id &id, const std::string& sceneName) const {
    return getItemProvider(getItemInfo(id, sceneName), sceneName);
}

// * * * * * * * * * * * * *  Public Member Functions  * * * * * * * * * * * //

ItemInfo* ItemStorage::createDeclareAndReturnCopyOfBaseItemInfo(const Id &id, const std::string& sceneName)
{
    const ItemInfo* baseInfo = getBaseItemInfo(id);
    ItemInfo* newInfo = baseInfo->copy();
    auto& itemInfos = sceneItemInfos_[sceneName].itemInfos;

    const auto& result = itemInfos.emplace(id, std::unique_ptr<ItemInfo>(newInfo));
    if (result.second)
    {
        return newInfo;
    }
    else
    {
        errorDeclaringItem(id, itemInfos, sceneName.c_str());
        return nullptr;
    }
}

bool ItemStorage::containsBaseItemInfo(const Id &id)
{
    return baseItemInfos_.itemInfos.find(id) != baseItemInfos_.itemInfos.end();
}

void ItemStorage::errorDeclaringItem(
    const Id &id,
    const std::unordered_map<Id, std::unique_ptr<ItemInfo>, IdHasher> &itemInfos,
    const char *sceneName
)
{
    std::string source = sceneName? (std::string("scene \"") + sceneName + "\"") : ("base");
    if (itemInfos.find(id) == itemInfos.end())
    {
        FatalErrorInFunction << "Unknown error trying to add item " << id.name << " to " << source.c_str()
                             << abort(FatalError);
    }
    else
    {
        FatalErrorInFunction << "Item " << id.name << " from " << source.c_str() << " was declared twice"
                             << abort(FatalError);
    }
}

const ItemInfo* ItemStorage::getBaseItemInfo(const Id &id) const
{
    const ItemInfo* result;
    try
    {
        result = baseItemInfos_.itemInfos.at(id).get();
    } catch (const std::exception& ex) {

        auto* modifiableBase = const_cast<SceneItemInfos*>(&baseItemInfos_);
        const auto& addResult = modifiableBase->itemInfos.emplace(id, std::unique_ptr<ItemInfo>(new MissingPostDictItemInfo(id)));
        result = addResult.first->second.get();
    }
    if (!result->isValid())
    {
        WarningInFunction << "Item " << id
                          << " was requested but " << result->getInvalidReason()
                          << endl;
    }
    return result;
}

const ItemInfo* ItemStorage::getItemInfo(const Id &id, const std::string& sceneName) const
{
    if (sceneItemInfos_.find(sceneName) != sceneItemInfos_.end())
    {
        auto iter = sceneItemInfos_.at(sceneName).itemInfos.find(id);
        if (iter != sceneItemInfos_.at(sceneName).itemInfos.end())
        {
            return iter->second.get();
        }
    }

    return getBaseItemInfo(id);
}

void ItemStorage::addToBoundariesListIfBoundary(const ItemInfo* itemInfo)
{
    const Id& id = itemInfo->getId();
    if (itemInfo->isMeshBoundary()) {
        meshBoundaries_(id.region).allBoundaries.push_back(id);
    }
    if (itemInfo->isExternalBoundary()) {
        meshBoundaries_(id.region).externalBoundaries.push_back(id);
    }
    if (itemInfo->isInternalBoundary()) {
        meshBoundaries_(id.region).internalBoundaries.push_back(id);
    }
    if (itemInfo->isProcessBoundary())
    {
        meshBoundaries_(id.region).processBoundaries.push_back(id);
    }
}

void ItemStorage::createMeshBoundaryGroup
(
        const std::string &name,
        const std::string &region,
        const std::vector<Id> &boundaryIds
)
{
    Id boundaryGroupId = Id(name, region, ItemType::Value::BOUNDARY_GROUP);
    if (baseItemInfos_.itemInfos.find(boundaryGroupId) == baseItemInfos_.itemInfos.end()) {
        if (!boundaryIds.empty()) {
            declareNewItemInfo(std::make_unique<BoundariesGroupInfo>(name, region, boundaryIds));
        }
    }
}

void ItemStorage::createMeshBoundaryGroups()
{
    for (const std::string& region : meshBoundaries_.toc())
    {
        createMeshBoundaryGroup(boundaryGroupKeys::EXTERNAL_BOUNDARY_GROUP_NAME, region, meshBoundaries_[region].externalBoundaries);
        createMeshBoundaryGroup(boundaryGroupKeys::INTERNAL_BOUNDARY_GROUP_NAME, region, meshBoundaries_[region].internalBoundaries);
    }
}

void ItemStorage::deleteProviders()
{
    auto it = providersFromItemInfo_.begin();
    while (it != providersFromItemInfo_.end())
    {
        if (it->first.itemInfo->isExternal())
        {
            it++;
        }
        else
        {
            it = providersFromItemInfo_.erase(it);
        }
    }
}

ItemInfo* ItemStorage::declareNewItemInfo(std::unique_ptr<ItemInfo> info)
{
    auto& itemInfos = baseItemInfos_.itemInfos;
    const Id id = info->getId();
    const auto& result = itemInfos.emplace(id, std::move(info));
    if (itemInfos.find(id) == itemInfos.end())
    {
        errorDeclaringItem(id, itemInfos);
    }
    addToBoundariesListIfBoundary(result.first->second.get());
    return result.first->second.get();
}

ItemInfo* ItemStorage::declareNewItemInfo(std::unique_ptr<ItemInfo> info, const std::string &sceneName)
{
    auto& itemInfos = sceneItemInfos_[sceneName].itemInfos;
    const Id id = info->getId();
    const auto& result = itemInfos.emplace(id, std::move(info));
    if (!result.second)
    {
        errorDeclaringItem(id, itemInfos, sceneName.c_str());
    }
    return result.first->second.get();
}


// * * * * * * * * * * * * *       Constructors        * * * * * * * * * * * //

ItemStorage::ItemStorage()
{
    sceneItemInfos_.clear();
    providersFromItemInfo_.clear();
}

// ************************************************************************* //
} // End namespace Foam
