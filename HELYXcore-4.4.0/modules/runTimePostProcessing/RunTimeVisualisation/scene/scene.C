/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.0.1
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
    (c) 2020-2023 Engys Ltd.

\*---------------------------------------------------------------------------*/

// OpenFOAM includes
#include "scene.H"

#include "db/Time/Time.H"
#include "rendering/compositers/compositerFactory.H"
#include "infos/runTimeVisualisationInfo.H"

#include "postDictObjectProviderDatabase.H"

// VTK includes
#include "vtkMultiProcessController.h"
#include "engysPVDWriter.h"

#include <utility>

namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

std::unique_ptr<Compositer> Scene::createCompositer
    (
        vtkMultiProcessController *controller,
        const RenderInfo &renderInfo,
        const HashTable<const ItemInfo *, Id, IdHasher> &itemInfos
    )
{
    bool hasTransparency = false;
    for (const ItemInfo *itemInfo: itemInfos)
    {
        if (itemInfo->isTransparent())
        {
            hasTransparency = true;
            break;
        }
    }
    return CompositerFactory::createCompositer(controller, hasTransparency, renderInfo);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Scene::Scene
    (
        std::string sceneName,
        const PostDictObjectProviderDatabase &database,
        const RunTimeVisualisationInfo &rtppInfo
    )
    :
    name_(std::move(sceneName)),
    sceneInfo_(database.getSceneInfo(name_)),
    database_(database),
    rtppInfo_(rtppInfo),
    compositer_(createCompositer(database.getController(), rtppInfo.getRenderInfo(), sceneInfo_.itemInfos)),
    lookupTableProvider_(sceneInfo_.colourLookupTablesInfo, database.getBaseColorLookupTablesInfo()),
    items_(database, sceneInfo_.itemInfos, name_),
    outputDirectory_(rtppInfo.getOutputFolder() / name_),
    renderManager_(
        sceneInfo_.backgroundColours,
        sceneInfo_.widgetsInfo,
        rtppInfo.getRenderInfo(),
        database_.getReferenceFrames(),
        compositer_.get()
    ),
    pvdWriter_(name_, rtppInfo, sceneInfo_, database.getController(), items_)
{
    items_.addToRenderer(renderManager_);
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Scene::renderAndExport
    (
        const Time &currentTime,
        const Time &runTime
    )
{
    if (!items_.isAnyVisible())
    {
        Warning << "Scene " << name_ << " has no visible items." << endl;
    }

    items_.updateActors(
        currentTime.timeIndex(), currentTime.value(),
        lookupTableProvider_,
        rtppInfo_.getColourMaps());

    pvdWriter_.exportIfRequired(currentTime);

    if (renderManager_.mustRenderImage() || renderManager_.mustRenderX3d())
    {
        if (renderManager_.mustRenderImage())
        {
            if (sceneInfo_.cameras.cameras.empty())
            {
                FatalError << "Scene " << name_ << " is set to export images, but no camera was defined" << endl << abort(FatalError);
            }
            int cameraIndex = 0;
            for (const CameraData& camera : sceneInfo_.cameras.cameras)
            {
                updateCamera(camera, currentTime, runTime);
                string cameraName = getCameraName(camera, cameraIndex, sceneInfo_.cameras.cameras.size());
                RendererExtraData extraData(
                    name_,
                    getOutputFileNameWithoutExtension(currentTime, cameraName),
                    currentTime.value(),
                    cameraName
                );
                renderManager_.renderAndExportPng(extraData);
                cameraIndex++;
            }
        }
        else
        {
            updateCamera(CameraData(), currentTime, runTime);
        }

        RendererExtraData extraData(
            name_,
            getOutputFileNameWithoutExtension(currentTime),
            currentTime.value(),
            ""
        );
        renderManager_.renderAndExportX3dIfNeeded(extraData);
        compositer_->resetForNewIteration();
    }
    Info<< decrIndent;
}

void Scene::updateCamera(const CameraData &camera, const Time &currentTime, const Time &runTime)
{
    renderManager_.setCamera(camera.toGlobalCamera());

    renderManager_.updateAndRedistributeWidgets(
        currentTime.value(),
        lookupTableProvider_,
        rtppInfo_.getColourMaps(),
        database_.getFoamMeshes(),
        database_.getExternalFields(),
        runTime
    );

    items_.registerDataToCompositer(compositer_.get());
    items_.redistributeDataWithCompositerIfNecessary(compositer_.get());
}

string Scene::getCameraName(const CameraData& camera, int cameraIndex, size_t totalCameras)
{
    if (totalCameras <= 1)
    {
        return "";
    }
    if (!camera.name.empty())
    {
        return camera.name;
    }
    return {"camera" + std::to_string(cameraIndex)};
}

string Scene::getOutputFileNameWithoutExtension(const Time &currentTime, const string& cameraName)
{
    if (cameraName.empty())
    {
        return getOutputFileNameWithoutExtension(currentTime);
    }
    else
    {
        return outputDirectory_ / name_ + "_" + cameraName + "_" + getTimeName(currentTime);
    }
}

string Scene::getOutputFileNameWithoutExtension(const Time &currentTime)
{
    return outputDirectory_ / name_ + "_" + getTimeName(currentTime);
}

string Scene::getTimeName(const Time &currentTime)
{
    // Ensure time steps have correct number of decimals
    const std::string::size_type requiredPrecision = currentTime.getTimeNamePrecision();
    word timeName = currentTime.timeName();
    std::string::size_type numberOfExtraChars = requiredPrecision - timeName.ext().size();

    if (timeName.ext().size() < requiredPrecision)
    {
        if (!timeName.hasExt())
        {
            timeName += ".";
        }
        timeName += string(numberOfExtraChars, '0');
    }
    else if (timeName.ext().size() > requiredPrecision)
    {
        WarningInFunction
            << "The user-defined time precision was set to "
            << requiredPrecision
            << ", which was insufficient to display the current time ("
            << timeName << ").  RunTimeVisualisation output files will be "
            << "numbered correctly, but due to the difference in number of "
            << "significant figures, your operating system may sort the "
            << "files in an unexpected order"
            << endl;
    }
    return timeName;
}


} // End namespace Foam

// ************************************************************************* //
