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
    Foam::functionObjects::runTimeVis::postDictObjectProviderDatabase

Group
    grpGraphicsFunctionObjects

Description
    Holds item data, which is shared between the function objects that use 3D objects.

    The functionality makes use of the VTK libraries (see http://www.vtk.org)
    which provide a broad set of functionality for scene composition and
    manipulation.

SourceFiles
    postDictObjectProviderDatabase.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "primitives/strings/word/word.H"
#include "db/Time/Time.H"
#include "baseClasses/foamField.H"
#include "infos/sceneInfos.H"
#include "storage/externalFields.H"

#include <memory>
#include <unordered_set>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class vtkMultiProcessController;

namespace Foam::functionObjects::runTimeVis
{
class Dictionaries;

/*---------------------------------------------------------------------------*\
                Class PostDictObjectProviderDatabase Declaration
\*---------------------------------------------------------------------------*/

struct PostDictObjectProviderDatabaseInstance;

class PostDictObjectProviderDatabase
{

private:

    static std::unique_ptr<PostDictObjectProviderDatabaseInstance> instance_;
    static unsigned int instanceReferenceCount_;
    static unsigned int staticReferenceCount_;
    static bool createdProviders_;
    static std::unordered_set<std::string> registeredFunctionObjects_;
    enum LoggerState
    {
        NOT_INITIALIZED,
        INITIALIZED,
        DEBUG
    };
    static LoggerState loggerInitialized_;
    static bool debug_;

    PostDictObjectProviderDatabaseInstance *localInstance_;

    word name;

public:

    PostDictObjectProviderDatabase(
        const word &name,
        const Time &runTime,
        const Dictionaries &dictionaries,
        bool debug
    );

    void addToRequiredItems(const Id &id);
    void addToRequiredScenes(const std::string &sceneName);

    void addToItemRequirements(const ItemRequirements &itemRequirements);

    size_t getSceneCount();
    [[nodiscard]] const SceneInfo &getSceneInfo(int index) const;
    [[nodiscard]] const SceneInfo &getSceneInfo(const std::string &sceneName) const;

    void updateMeshDomainForTimestep(label timeIndex, scalar currentTime);
    void updateExternalDomainForTimestep(label timeIndex);

    [[nodiscard]] vtkSmartPointer<vtkDataSet> getDataSetForSceneItem(
        const ItemInfo *info,
        const std::string &sceneName,
        label timeIndex,
        scalar currentTime
    ) const;
    [[nodiscard]] vtkSmartPointer<vtkPolyData> getOutlineForSceneItem(
        const ItemInfo *info,
        const std::string &sceneName,
        label timeIndex,
        scalar currentTime
    ) const;
    [[nodiscard]] vtkSmartPointer<vtkDataSet> getDataSetForBaseItem(
        const Id &id,
        label timeIndex,
        scalar currentTime
    ) const;

    [[nodiscard]] const ColourLookupTablesInfo &getBaseColorLookupTablesInfo() const;

    [[nodiscard]] const FoamMeshes &getFoamMeshes() const;
    [[nodiscard]] const ExternalFields &getExternalFields() const;

    [[nodiscard]] vtkMultiProcessController *getController() const;

    void createItemProviders(const Time& time);
    void deleteItemProviders(const Time& time);

    [[nodiscard]] const ReferenceFrames& getReferenceFrames() const;

    static void initializeStaticCommon(bool debug);

    static void finalizeStaticCommon();

    ~PostDictObjectProviderDatabase();

private:
    static void initialiseMpi();
    static void initializeLogger(bool debug);
    static bool isLastFunctionObject(const word& name, const Time &time);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
