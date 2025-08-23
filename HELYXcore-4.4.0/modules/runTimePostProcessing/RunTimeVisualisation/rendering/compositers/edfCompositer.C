/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : Dev
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

#include "compileOptions.H"

// OpenFOAM includes
#include "edfCompositer.H"

#include "rendering/rtppActor.H"
#include "rendering/edfRenderers/edfWriter.H"
#include "infos/renderInfo.H"

// vtk includes
#include "vtkCompositedSynchronizedRenderers.h"
#include "vtkRenderWindow.h"
#include "engysEDFOpenGLRenderer.h"
#include "engysEDFPropsToRender.h"
#include "engysEDFField.h"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

    EDFCompositer::EDFCompositer(vtkMultiProcessController* controller, const RenderInfo& renderInfo)
:
    Compositer(controller)
{
    vtkNew<engysEDFOpenGLRenderer> renderer;
    renderer->SetController(controller);
    this->propsToRender_ = vtkSmartPointer<engysEDFPropsToRender>::New();
    this->useFieldsShownByProps = renderInfo.edfFields.empty();
    for (const foamField &field: renderInfo.edfFields)
    {
        vtkNew<engysEDFField> engysField;
        if (engysField->FromFoamField(field.c_str()))
        {
            this->propsToRender_->AddFieldToAllProps(engysField);
        }
    }
    renderer->SetPropsToRender(propsToRender_);
    setRendererWriter(new EdfWriter(renderInfo, renderer));
    renderer_ = renderer;
    synchronizer_ = vtkSmartPointer<vtkSynchronizedRenderers>::New();
    synchronizer_->ParallelRenderingOn();
    synchronizer_->SetParallelController(controller_);
}

void EDFCompositer::registerActorPolyData(rtppActor& actor, const char* name)
{
    std::vector<vtkActor*> vtkActors = actor.getVTKActors();
    std::vector<std::string> vtkActorSuffixes = actor.getVTKActorNameSuffixes();
    for (size_t i = 0; i < vtkActors.size(); i++)
    {
        std::string vtkActorName = name + vtkActorSuffixes[i];
        this->propsToRender_->AddProp(vtkActors[i], vtkActorName.c_str());
        if (this->useFieldsShownByProps)
        {
            this->propsToRender_->AddFieldShownByPropToProp(vtkActors[i]);
        }
    }
}

void EDFCompositer::renderAndWriteImage
(
        vtkRenderWindow* window,
        const RendererExtraData& extraData
)
{
    this->synchronizer_->SetRenderer(this->renderer_);
    vtkSmartPointer<vtkSynchronizedRenderWindows> windowSync = createWindowSynchronizer(window);
    renderAndWriteImage_(window, extraData);
}

void EDFCompositer::resetForNewIteration()
{
    this->propsToRender_->ClearProps();
}

} // End namespace Foam
