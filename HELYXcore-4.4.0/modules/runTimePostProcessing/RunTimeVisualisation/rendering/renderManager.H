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
    (c) 2020-2021 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::renderManager

Description
    Class to control rendering

SourceFiles
    renderManager.C

\*---------------------------------------------------------------------------*/

#pragma once

// OpenFOAM includes
#include "primitives/ints/label/label.H"
#include "infos/sceneInfo.H"
#include "infos/renderInfo.H"
#include "widgets/widgetsInScene.H"
#include "widgets/widgetsInOverlay.H"
#include "widgets/overlay/colorLegend.H"
#include "widgets/overlay/axesWidget.H"
#include "widgets/overlay/logoWidget.H"
#include "widgets/overlay/timestepWidget.H"
#include "rendering/rendererExtraData.H"

// VTK includes
#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward declarations
class vtkRenderer;
class vtkRenderWindow;
class vtkRenderWindowInteractor;
class vtkMultiProcessController;
class vtkActor;
class vtkOpenGLPolyDataMapper;

namespace Foam::functionObjects::runTimeVis
{

class ColourLookupTableProvider;
class Compositer;
/*---------------------------------------------------------------------------*\
                           Class renderManager Declaration
\*---------------------------------------------------------------------------*/

class renderManager
{
    // Private data

        vtkSmartPointer<vtkRenderer> baseRenderer_;
        vtkSmartPointer<vtkRenderer> overlayRenderer_;
        vtkSmartPointer<vtkRenderWindow> renderWindow_;
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor_;
        std::vector<vtkOpenGLPolyDataMapper*> mappers_;

        RenderInfo renderInfo_;

        WidgetsInScene widgetsInScene_;
        WidgetsInOverlay widgetsInOverlay_;

        Compositer* compositer_;

    // Private Member Functions

        void initialiseRenderer(const BackgroundColourData& background);

public:

    // Constructors

        renderManager
        (
                const BackgroundColourData& background,
                const WidgetsInfo& widgetsInfo,
                const RenderInfo& renderInfo,
                const ReferenceFrames& referenceFrames,
                Compositer* compositer
        );

        void addActor(rtppActor& actor);


    // Member Functions
        bool mustRenderImage() const;
        void renderAndExportPng(const RendererExtraData& extraData);
        bool mustRenderX3d() const;
        void renderAndExportX3dIfNeeded(const RendererExtraData& extraData);

        void setCamera(const CameraData& cameraData);

        void updateAndRedistributeWidgets
        (
                scalar currentTimeValue,
                ColourLookupTableProvider& colourLutProvider,
                const ColourMaps& colourMaps,
                const FoamMeshes& meshes,
                const ExternalFields& externalFields,
                const Time& runTime
        );


    // Access

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
