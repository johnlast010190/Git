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
    (c) 2015-2019 OpenCFD Ltd.
    (c) 2020-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "postDictObjectProviderDatabaseInstance.H"

#include "infos/items/baseItemInfoReader.H"
#include "Utils/ParallelUtils.H"
#include "dictionaries/dictionaries.H"

#include "primitives/strings/fileName/fileName.H"

static const char *const BASE_SCENE_NAME = "__baseModel__";

namespace Foam::functionObjects::runTimeVis
{

PostDictObjectProviderDatabaseInstance::PostDictObjectProviderDatabaseInstance
    (
        const Time &runTime,
        const Dictionaries &dictionaries
    ) :
    itemStorage_(),
    meshes_(runTime, dictionaries),
    externalFields_(),
    referenceFrames_(dictionaries, meshes_.getDefaultMesh().getMesh()),
    factory_((ParallelUtils::isRunningInParallel() ? runTime.path() / ".." : runTime.path()), meshes_, itemStorage_, referenceFrames_, dictionaries),
    baseColorLookupTablesInfo_(dictionaries.getRtppDict()),
    instanceServices_(*this)
{
    updateMeshDomainForTimestep(-1, -1);

    BaseItemInfoReader::readBaseItemsToFactory(factory_, meshes_, dictionaries);

    sceneInfos_.readDictionaryAndCreateItemInfos(dictionaries.getRtppDict(), factory_, referenceFrames_);
    addFieldsToMeshFieldList(dictionaries.getRtppDict(), sceneInfos_);
}

void PostDictObjectProviderDatabaseInstance::addFieldsToMeshFieldList(const dictionary &postDict, const SceneInfos& sceneInfos)
{
    const dictionary &legendsDict = GET_OPTIONAL_DICTIONARY(postDict, sceneKeys::COLOR_LEGENDS_DICT_KEY);

    for (const std::string &fieldName: legendsDict.toc())
    {
        meshes_.addColourField(foamField(fieldName));
        externalFields_.addColourField(foamField(fieldName));
    }
    for (const foamField &fieldName: sceneInfos.getItemRequirements().getRequiredFields().getColorFields())
    {
        meshes_.addColourField(foamField(fieldName));
        externalFields_.addColourField(foamField(fieldName));
    }
}

void PostDictObjectProviderDatabaseInstance::updateMeshDomainForTimestep(label timeIndex, scalar timeValue)
{
    instanceServices_.updateIfNecessary(timeIndex, timeValue);
    if (timeIndex > timeIndexOfLastMeshUpdate)
    {
        meshes_.updateDomainBounds();
        meshes_.updateDomainRangeForColourFields();
        timeIndexOfLastMeshUpdate = timeIndex;
    }
}

void PostDictObjectProviderDatabaseInstance::updateExternalDomainForTimestep(label timeIndex)
{
    if (timeIndex > timeIndexOfLastExternalUpdate)
    {
        externalFields_.updateForTimeStep(timeIndex);
        timeIndexOfLastExternalUpdate = timeIndex;
    }
}

void PostDictObjectProviderDatabaseInstance::addToRequiredItems(const Id &id)
{
    baseRequiredItems_.push_back(itemStorage_.getBaseItemInfo(id));
}

void PostDictObjectProviderDatabaseInstance::addToRequiredItems(const ItemInfo &itemInfo, const std::string &scene)
{
    sceneRequiredItems_[scene].push_back(&itemInfo);
}

void PostDictObjectProviderDatabaseInstance::addToRequiredScenes(const std::string &sceneName)
{
    const SceneInfo &sceneInfo = getSceneInfo(sceneName);
    std::for_each(
        sceneInfo.itemInfos.begin(), sceneInfo.itemInfos.end(),
        [this, sceneName](const ItemInfo *itemInfo) {
            this->addToRequiredItems(*itemInfo, sceneName);
        }
    );
}

void PostDictObjectProviderDatabaseInstance::addToItemRequirements(const ItemRequirements &itemRequirements)
{
    extraRequirements_.mergeFromDownstream(itemRequirements);
}

size_t PostDictObjectProviderDatabaseInstance::getSceneCount() const
{
    return sceneInfos_.getSceneCount();
}

const SceneInfo &PostDictObjectProviderDatabaseInstance::getSceneInfo(int index) const
{
    return sceneInfos_.getSceneInfo(index);
}

const SceneInfo &PostDictObjectProviderDatabaseInstance::getSceneInfo(const std::string &sceneName) const
{
    return sceneInfos_.getSceneInfo(sceneName);
}

vtkMultiProcessController *PostDictObjectProviderDatabaseInstance::getController() const
{
    return vtkMultiProcessController::GetGlobalController();
}

ItemDataSetProvider *PostDictObjectProviderDatabaseInstance::getDataSetProviderForSceneItem(
    const ItemInfo *info,
    const std::string &sceneName
) const
{
    return itemStorage_.getItemProvider(info, sceneName);
}

ItemDataSetProvider *PostDictObjectProviderDatabaseInstance::getDataSetProviderForBaseItem(const Id &id) const
{
    return itemStorage_.getItemProvider(id, BASE_SCENE_NAME);
}

const ColourLookupTablesInfo &PostDictObjectProviderDatabaseInstance::getBaseColorLookupTablesInfo() const
{
    return baseColorLookupTablesInfo_;
}

const FoamMeshes &PostDictObjectProviderDatabaseInstance::getFoamMeshes() const
{
    return meshes_;
}

void PostDictObjectProviderDatabaseInstance::createItemProviders()
{
    createItemProviders(BASE_SCENE_NAME, baseRequiredItems_);

    std::for_each(
        sceneRequiredItems_.begin(), sceneRequiredItems_.end(),
        [this](const std::pair<std::string, std::vector<const ItemInfo *>> &scenePair) {
            createItemProviders(scenePair.first, scenePair.second);
        }
    );
}

void PostDictObjectProviderDatabaseInstance::createItemProviders(
    const std::string &sceneName,
    const std::vector<const ItemInfo *> &sceneItems
)
{
    std::for_each(
        sceneItems.begin(), sceneItems.end(),
        [this, sceneName](const ItemInfo *item) {
            ItemDataSetProvider *provider = itemStorage_.createItemProvider(
                item,
                sceneName,
                extraRequirements_,
                &instanceServices_
            );
            if (item->isExternal())
            {
                externalFields_.insertProvider(static_cast<ExternalItemDataSetProvider*>(provider));
            }
        }
    );
}

void PostDictObjectProviderDatabaseInstance::deleteItemProviders()
{
    itemStorage_.deleteProviders();
    externalFields_.clearProviders();
}

PostDictObjectProviderDatabaseInstance::~PostDictObjectProviderDatabaseInstance()
{
    itemStorage_.deleteProviders();
}

}

// ************************************************************************* //
