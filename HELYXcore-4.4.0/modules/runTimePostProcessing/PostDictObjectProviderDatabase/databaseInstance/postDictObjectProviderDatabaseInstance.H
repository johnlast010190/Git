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

Class
    Foam::functionObjects::runTimeVis::postDictObjectProviderDatabaseInstance

Group
    grpGraphicsFunctionObjects

Description
    Holds item data, which is shared between the function objects that use 3D objects.

    The functionality makes use of the VTK libraries (see http://www.vtk.org)
    which provide a broad set of functionality for scene composition and
    manipulation.

SourceFiles
    postDictObjectProviderDatabaseInstance.C

\*---------------------------------------------------------------------------*/

#pragma once

#include <memory>

#include "infos/sceneInfos.H"
#include "storage/foamMeshes.H"
#include "storage/externalFields.H"
#include "storage/itemStorage.H"
#include "storage/referenceFrames.H"

#include "instanceServices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward declarations
class vtkMultiProcessController;

namespace Foam::functionObjects::runTimeVis
{
class Dictionaries;

/*---------------------------------------------------------------------------*\
          Class postDictObjectProviderDatabaseInstance Declaration
\*---------------------------------------------------------------------------*/

class PostDictObjectProviderDatabaseInstance
{
private:

    ItemStorage itemStorage_;
    FoamMeshes meshes_;
    ExternalFields externalFields_;
    ReferenceFrames referenceFrames_;
    ItemInfoFactory factory_;
    SceneInfos sceneInfos_;

    const ColourLookupTablesInfo baseColorLookupTablesInfo_;

    std::unordered_map<std::string, std::vector<const ItemInfo *>> sceneRequiredItems_;
    std::vector<const ItemInfo *> baseRequiredItems_;

    InstanceServices instanceServices_;

    ItemRequirements extraRequirements_;

    label timeIndexOfLastMeshUpdate = labelMin;
    label timeIndexOfLastExternalUpdate = labelMin;

public:
    PostDictObjectProviderDatabaseInstance(
        const Time &runTime,
        const Dictionaries &dictionaries
    );

    void addToRequiredItems(const Id &id);
    void addToRequiredItems(const ItemInfo &itemInfo, const std::string &scene);
    void addToRequiredScenes(const std::string &sceneName);

    void addToItemRequirements(const ItemRequirements &itemRequirements);

    [[nodiscard]] size_t getSceneCount() const;
    [[nodiscard]] const SceneInfo &getSceneInfo(int index) const;
    [[nodiscard]] const SceneInfo &getSceneInfo(const std::string &sceneName) const;

    [[nodiscard]] vtkMultiProcessController *getController() const;

    [[nodiscard]] ItemDataSetProvider *getDataSetProviderForSceneItem(const ItemInfo *info, const std::string &sceneName) const;
    [[nodiscard]] ItemDataSetProvider *getDataSetProviderForBaseItem(const Id &id) const;

    [[nodiscard]] const ColourLookupTablesInfo &getBaseColorLookupTablesInfo() const;

    [[nodiscard]] const FoamMeshes &getFoamMeshes() const;
    [[nodiscard]] const ExternalFields &getExternalFields() const
    {
        return externalFields_;
    }

    void updateMeshDomainForTimestep(label timeIndex, scalar timeValue);
    void updateExternalDomainForTimestep(label timeIndex);

    void createItemProviders();
    void deleteItemProviders();

    [[nodiscard]] const ReferenceFrames& getReferenceFrames() const
    {
        return referenceFrames_;
    }

    ~PostDictObjectProviderDatabaseInstance();

private:
    void addFieldsToMeshFieldList(const dictionary &postDict, const SceneInfos& sceneInfos);
    void createItemProviders(
        const std::string &sceneName,
        const std::vector<const ItemInfo *> &sceneItems
    );

    friend class InstanceServicesRegion;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
