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

Class
    Foam::functionObjects::runTimeVis::ItemStorage

Description
    Stores all the items read from the dictionaries, allowing them to be accessed
    with a combination of their ID and scene

SourceFiles
    itemStorage.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "infos/items/itemInfo.H"
#include "storage/itemInfoIndex.H"
#include "baseClasses/itemRequirements.H"

#include <unordered_map>
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{
class InstanceServices;

/*---------------------------------------------------------------------------*\
                        Class ItemStorage Declaration
\*---------------------------------------------------------------------------*/

class ItemStorage
{
private:

    struct SceneItemInfos
    {
        std::unordered_map<Id, std::unique_ptr<ItemInfo>, IdHasher> itemInfos;

        SceneItemInfos()
        {
            itemInfos.clear();
        };
    };

    struct MeshBoundaryIds
    {
        std::vector<Id> allBoundaries;
        std::vector<Id> externalBoundaries;
        std::vector<Id> internalBoundaries;
        std::vector<Id> processBoundaries;
    };

    // Private data

    // The item storage needs the infos from all scenes so that identical items can be merged
    std::unordered_map<std::string, SceneItemInfos> sceneItemInfos_;
    SceneItemInfos baseItemInfos_;

    std::unordered_map<ItemInfoIndex, std::unique_ptr<ItemDataSetProvider>, ItemInfoIndexHasher> providersFromItemInfo_;

    HashTable<MeshBoundaryIds, std::string> meshBoundaries_;

    // Private Member Functions

    void addSourceToProvider(ItemDataSetProvider *provider, const Id &sourceId, const std::string &sceneName);
    void addSourcesToProvider(ItemDataSetProvider *provider, const ItemInfo *itemInfo, const std::string &sceneName);

    void createMeshBoundaryGroup
        (
            const std::string &name,
            const std::string &region,
            const std::vector<Id> &boundaryIds
        );

    void addToBoundariesListIfBoundary(const ItemInfo *info);

    [[nodiscard]] ItemDataSetProvider *getOrCreateItemProvider(const ItemInfo *itemInfo, const std::string &sceneName);
    [[nodiscard]] ItemDataSetProvider *getItemProviderIfExists(const ItemInfoIndex &itemInfoIndex) const;

public:

    // Constructors

    //- Create an empty storage
    ItemStorage();

    // Member Functions

    ItemInfo* declareNewItemInfo(std::unique_ptr<ItemInfo> info);
    ItemInfo* declareNewItemInfo(std::unique_ptr<ItemInfo> info, const std::string &sceneName);

    ItemInfo *createDeclareAndReturnCopyOfBaseItemInfo(const Id &id, const std::string &sceneName);

    [[nodiscard]] const ItemInfo *getBaseItemInfo(const Id &id) const;

    [[nodiscard]] const ItemInfo *getItemInfo(const Id &id, const std::string &sceneName) const;

    ItemDataSetProvider* createItemProvider(
        const ItemInfo *itemInfo,
        const std::string &sceneName,
        const ItemRequirements& extraRequirements,
        InstanceServices* instanceServices
    );

    ItemDataSetProvider *getItemProvider(const ItemInfo *itemInfo, const std::string &sceneName) const;
    ItemDataSetProvider *getItemProvider(const Id &id, const std::string& sceneName) const;

    void createMeshBoundaryGroups();

    std::vector<ItemDataSetProvider*> getAllMeshBoundaryProviders(const std::string& region);
    std::vector<ItemDataSetProvider*> getProcessBoundaryProviders(const std::string& region);

    void deleteProviders();

    [[nodiscard]] bool containsBaseItemInfo(const Id &id);

private:
    static void errorDeclaringItem(
        const Id &id,
        const std::unordered_map<Id, std::unique_ptr<ItemInfo>, IdHasher> &itemInfos,
        const char *sceneName = nullptr
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
