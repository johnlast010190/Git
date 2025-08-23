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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "axesWidget.H"

#include "dataStructs/widgets/axisWidgetData.H"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkOrientationMarkerWidget.h"
#include "vtkAxesActor.h"
#include "vtkCaptionActor2D.h"
#include "vtkTextProperty.h"
#include "vtkTextActor.h"
#include "vtkRenderWindowInteractor.h"

namespace Foam::functionObjects::runTimeVis
{

std::unique_ptr<AxesWidget> AxesWidget::createWidgetIfVisible(const AxisWidgetData& data)
{
    if (data.visible)
    {
        return std::unique_ptr<AxesWidget>(new AxesWidget(data));
    }
    else
    {
        return {nullptr};
    }
}


AxesWidget::AxesWidget(const AxisWidgetData&)
{
    auto axes = vtkSmartPointer<vtkAxesActor>::New();
    axes->GetXAxisCaptionActor2D()->GetTextActor()->GetTextProperty()->ShadowOff();
    axes->GetXAxisCaptionActor2D()->GetTextActor()->GetTextProperty()->SetFontFamilyToArial();
    axes->GetXAxisCaptionActor2D()->GetTextActor()->GetTextProperty()->ItalicOff();
    axes->GetXAxisCaptionActor2D()->GetTextActor()->GetTextProperty()->BoldOff();
    axes->GetXAxisCaptionActor2D()->GetTextActor()->GetTextProperty()->SetColor(0, 0, 0);

    axes->GetYAxisCaptionActor2D()->GetTextActor()->GetTextProperty()->ShadowOff();
    axes->GetYAxisCaptionActor2D()->GetTextActor()->GetTextProperty()->SetFontFamilyToArial();
    axes->GetYAxisCaptionActor2D()->GetTextActor()->GetTextProperty()->ItalicOff();
    axes->GetYAxisCaptionActor2D()->GetTextActor()->GetTextProperty()->BoldOff();
    axes->GetYAxisCaptionActor2D()->GetTextActor()->GetTextProperty()->SetColor(0, 0, 0);

    axes->GetZAxisCaptionActor2D()->GetTextActor()->GetTextProperty()->ShadowOff();
    axes->GetZAxisCaptionActor2D()->GetTextActor()->GetTextProperty()->SetFontFamilyToArial();
    axes->GetZAxisCaptionActor2D()->GetTextActor()->GetTextProperty()->ItalicOff();
    axes->GetZAxisCaptionActor2D()->GetTextActor()->GetTextProperty()->BoldOff();
    axes->GetZAxisCaptionActor2D()->GetTextActor()->GetTextProperty()->SetColor(0, 0, 0);

    axesWidget_ = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
    axesWidget_->SetOrientationMarker(axes);
    axesWidget_->SetViewport(0, 0, 0.25, 0.25);
}

void AxesWidget::addToRenderer(vtkRenderer* renderer)
{
    if (!axesWidget_->GetInteractor())
    {
        axesWidget_->SetInteractor(renderer->GetRenderWindow()->GetInteractor());
        axesWidget_->SetCurrentRenderer(renderer);
    }
    axesWidget_->On();
}

void AxesWidget::removeFromRenderer()
{
    axesWidget_->Off();
}

// ************************************************************************* //
} // End namespace
