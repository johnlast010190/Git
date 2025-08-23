/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.3.0
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
    (c) 2020-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

// OpenFOAM includes
#include "renderManager.H"
#include "Utils/ParallelUtils.H"
#include "Utils/boundsUtils.H"
#include "storage/foamMeshes.H"
#include "colourLookupTable/colourLookupTableProvider.H"
#include "rendering/compositers/compositer.H"
#include "primitives/Scalar/scalar/scalar.H"
#include "rendering/rtppActor.H"

// VTK includes
#include "vtkLight.h"
#include "vtkCamera.h"
#include "vtkProperty2D.h"
#include "vtkActor.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkProperty.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkOpenGLPolyDataMapper.h"

enum Layers {
    BASE_LAYER,
    OVERLAY_LAYER,
    NUMBER_OF_LAYERS
};

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

vtkSmartPointer<vtkLight> createLight(Foam::scalar elevation, Foam::scalar azimuth)
{
    auto light = vtkSmartPointer<vtkLight>::New();
    light->SetIntensity(0.5);
    light->SetColor(1, 1, 1);
    light->SetLightTypeToCameraLight();
    light->SetDirectionAngle(elevation, azimuth);
    light->SwitchOn();

    return light;
}

namespace Foam::functionObjects::runTimeVis
{

void renderManager::setCamera(const CameraData &cameraData)
{
    vtkNew<vtkCamera> camera;
    camera->SetParallelProjection
        (
            cameraData.parallelProjection
        );


    if (cameraData.parallelProjection)
    {
        camera->SetParallelScale(cameraData.parallelScale);
    }

    camera->SetViewUp
        (
            cameraData.up.x(),
            cameraData.up.y(),
            cameraData.up.z()
        );
    camera->SetPosition
        (
            cameraData.position.x(),
            cameraData.position.y(),
            cameraData.position.z()
        );
    camera->SetFocalPoint
        (
            cameraData.focalPoint.x(),
            cameraData.focalPoint.y(),
            cameraData.focalPoint.z()
        );

    camera->Modified();
    baseRenderer_->SetActiveCamera(camera);
    overlayRenderer_->SetActiveCamera(camera);
    // Computing normals in the cpu instead of on the shader is only faster if we are
    // rendering several frames, which is not the case, and only gives better results
    // if the camera is parallel
    for (vtkOpenGLPolyDataMapper *mapper: mappers_)
    {
        mapper->SetComputeNormalsInCPU(cameraData.parallelProjection);
    }
}


void renderManager::initialiseRenderer(const BackgroundColourData &background)
{
    renderWindow_ = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow_->SetNumberOfLayers(NUMBER_OF_LAYERS);
    renderWindow_->OffScreenRenderingOn();
    renderWindow_->SetMultiSamples(0);
    if (renderInfo_.exportPng || renderInfo_.exportEdf)
    {
        renderWindow_->SetSize(renderInfo_.width, renderInfo_.height);
    }

    renderWindowInteractor_ = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor_->SetRenderWindow(renderWindow_);

    baseRenderer_->SetLayer(BASE_LAYER);

    baseRenderer_->AddLight(createLight(45, 45));
    baseRenderer_->AddLight(createLight(-45, 45));
    baseRenderer_->AddLight(createLight(45, -45));
    baseRenderer_->AddLight(createLight(-45, -45));

    baseRenderer_->SetUseDepthPeeling(true);
    baseRenderer_->SetMaximumNumberOfPeels(4);
    baseRenderer_->SetOcclusionRatio(0);

    // postDict always sets both background 1 and 2, but uses them
    // backwards...
    if (renderInfo_.transparentBackground)
    {
        baseRenderer_->GradientBackgroundOff();
        baseRenderer_->SetBackgroundAlpha(0);
    }
    else
    {
        baseRenderer_->GradientBackgroundOn();
        baseRenderer_->SetBackgroundAlpha(1);
        baseRenderer_->SetBackground
            (
                background.background2.x(),
                background.background2.y(),
                background.background2.z()
            );
        baseRenderer_->SetBackground2
            (
                background.background1.x(),
                background.background1.y(),
                background.background1.z()
            );
    }
    renderWindow_->AddRenderer(baseRenderer_);

    overlayRenderer_ = vtkSmartPointer<vtkRenderer>::New();
    overlayRenderer_->SetLayer(OVERLAY_LAYER);
    renderWindow_->AddRenderer(overlayRenderer_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

renderManager::renderManager
    (
        const BackgroundColourData &background,
        const WidgetsInfo &widgetsInfo,
        const RenderInfo &renderInfo,
        const ReferenceFrames& referenceFrames,
        Compositer *compositer
    )
    :
    baseRenderer_(compositer->get3DRenderer()),
    renderInfo_(renderInfo),
    widgetsInScene_(compositer, widgetsInfo, referenceFrames),
    widgetsInOverlay_(widgetsInfo),
    compositer_(compositer)
{
    initialiseRenderer(background);
    mappers_.clear();

    vtkMapper::SetResolveCoincidentTopologyToPolygonOffset();
    vtkMapper::SetResolveCoincidentTopologyPolygonOffsetParameters(0, 0);
    vtkMapper::SetResolveCoincidentTopologyLineOffsetParameters(-0.5, 0);
    vtkMapper::SetResolveCoincidentTopologyPointOffsetParameter(-2);

    widgetsInScene_.setRenderer(baseRenderer_);
    widgetsInScene_.addVisibleWidgetsToRenderer();

    widgetsInOverlay_.setRenderer(overlayRenderer_);
    widgetsInOverlay_.addVisibleWidgetsToRenderer();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::runTimeVis::renderManager::addActor(rtppActor& rtppActor)
{
    std::vector<vtkActor*> vtkActors = rtppActor.getVTKActors();
    for (vtkActor* actor : vtkActors)
    {
        if (actor)
        {
            mappers_.push_back(vtkOpenGLPolyDataMapper::SafeDownCast(actor->GetMapper()));
            baseRenderer_->AddActor(actor);
        }
    }
}

void renderManager::updateAndRedistributeWidgets
    (
        scalar currentTimeValue,
        ColourLookupTableProvider &colourLutProvider,
        const ColourMaps &colourMaps,
        const FoamMeshes &meshes,
        const ExternalFields& externalFields,
        const Time &runTime
    )
{
    widgetsInOverlay_.updateVisibleWidgets(currentTimeValue, colourLutProvider, colourMaps, meshes, externalFields, runTime);

    widgetsInScene_.updateVisibleWidgets();
    widgetsInScene_.registerDataToCompositer();
    widgetsInScene_.redistributeDataWithCompositer();
}

bool renderManager::mustRenderImage() const
{
    return renderInfo_.exportPng || renderInfo_.exportEdf;
}

void Foam::functionObjects::runTimeVis::renderManager::renderAndExportPng(const RendererExtraData& extraData)
{
    if (renderInfo_.exportPng || renderInfo_.exportEdf)
    {
        boundBox bounds = boundsUtils::computeAllProcsRenderBoundingBox(baseRenderer_);
        scalar sBounds[6];
        boundsUtils::scalarArrayFromBoundBox(sBounds, bounds);
        baseRenderer_->ResetCameraClippingRange(sBounds[0], sBounds[1], sBounds[2], sBounds[3], sBounds[4], sBounds[5]);

        compositer_->renderAndWriteImage(renderWindow_, extraData);
    }
    if (renderInfo_.exportX3d)
    {
        compositer_->renderAndWriteX3D(renderWindow_, extraData.pathToOutputFileWithoutExtension);
    }
}

bool renderManager::mustRenderX3d() const
{
    return renderInfo_.exportX3d;
}

void renderManager::renderAndExportX3dIfNeeded(const RendererExtraData& extraData)
{
    if (mustRenderX3d())
    {
        compositer_->renderAndWriteX3D(renderWindow_, extraData.pathToOutputFileWithoutExtension);
    }
}

} // namespace
// ************************************************************************* //
