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

#include "gridWidget.H"

#include "Utils/boundsUtils.H"

#include <utility>

#include "vtkGridAxes3DActor.h"
#include "vtkTextProperty.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkProp.h"
#include "vtkMapper.h"

namespace Foam::functionObjects::runTimeVis
{

GridWidget::GridWidget(GridWidgetData  data)
:
        data_(std::move(data))
{
    if (!isVisible()) {
        return;
    }

    axesActor_ = vtkSmartPointer<vtkGridAxes3DActor>::New();

    configureGridColor();
    configureLabelColor();
    configureAxesTitles();
    configureFaces();
    configureFonts();
    configureCustomBounds();
}

#define PassColor(color) (color)[0], (color)[1], (color)[2]

void GridWidget::configureGridColor() {
    vtkNew<vtkProperty> property;
    property->DeepCopy(axesActor_->GetProperty());
    property->SetColor(PassColor(data_.axesGridColor));
    axesActor_->SetProperty(property);
}

void GridWidget::configureLabelColor() {
    for(int axis = 0; axis < 3; axis++) {
        axesActor_->GetTitleTextProperty(axis)->SetColor(PassColor(data_.axesGridLabelColor));
        axesActor_->GetLabelTextProperty(axis)->SetColor(PassColor(data_.axesGridLabelColor));
    }
}

void GridWidget::configureAxesTitles() {
    if (data_.showAxesGridTitle)
    {
        axesActor_->SetTitle(0, data_.XAxisTitle);
        axesActor_->SetTitle(1, data_.YAxisTitle);
        axesActor_->SetTitle(2, data_.ZAxisTitle);
    }
    else
    {
        axesActor_->SetTitle(0, "");
        axesActor_->SetTitle(1, "");
        axesActor_->SetTitle(2, "");
    }
}

#define ALL_FACES_MASK UINT_MAX

void GridWidget::configureFaces() {
    vtkNew<vtkProperty> property;
    property->DeepCopy(axesActor_->GetProperty());
    if (data_.manualFaceSelect)
    {
        axesActor_->SetFaceMask(getFaceMaskFromData());
        property->SetFrontfaceCulling(false);
    }
    else
    {
        axesActor_->SetFaceMask(ALL_FACES_MASK);
        property->SetFrontfaceCulling(true);
    }
    axesActor_->SetProperty(property);
}

unsigned int GridWidget::getFaceMaskFromData() const {
    return ((data_.showXYMin ? vtkGridAxes3DActor::MIN_XY : 0) |
            (data_.showYZMin ? vtkGridAxes3DActor::MIN_YZ : 0) |
            (data_.showZXMin ? vtkGridAxes3DActor::MIN_ZX : 0) |
            (data_.showXYMax ? vtkGridAxes3DActor::MAX_XY : 0) |
            (data_.showYZMax ? vtkGridAxes3DActor::MAX_YZ : 0) |
            (data_.showZXMax ? vtkGridAxes3DActor::MAX_ZX : 0));
}

void GridWidget::configureFonts() {
    for (int axis = 0; axis < 3; axis++) {
        axesActor_->GetTitleTextProperty(axis)->SetFontFamilyToArial();
        axesActor_->GetTitleTextProperty(axis)->SetFontSize(FONT_SIZE);
        axesActor_->GetLabelTextProperty(axis)->SetFontFamilyToArial();
        axesActor_->GetLabelTextProperty(axis)->SetFontSize(FONT_SIZE);
    }
}

void GridWidget::configureCustomBounds() {
    if (data_.useCustomBounds) {
        scalar bounds[6] = { data_.XBounds[0], data_.XBounds[1],
                data_.YBounds[0], data_.YBounds[1], data_.ZBounds[0],
                data_.ZBounds[1] };
        configureCustomLabels(bounds);
        axesActor_->SetGridBounds(bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);
    }
}

void GridWidget::configureCustomLabels(scalar bounds[6]) {
    if (data_.useCustomLabels) {
        configureLabelForAxis(bounds, data_.useXDistanceLabel, data_.XDistanceLabel, X);
        configureLabelForAxis(bounds, data_.useYDistanceLabel, data_.YDistanceLabel, Y);
        configureLabelForAxis(bounds, data_.useZDistanceLabel, data_.ZDistanceLabel, Z);
    } else {
        removeLabels();
    }
}

void GridWidget::configureLabelForAxis
(
        const scalar bounds[6],
        bool useLabelDistance,
        scalar labelDistance,
        label axis
)
{
    if (useLabelDistance && labelDistance > 0)
    {
        axesActor_->SetUseCustomLabels(axis, true);

        scalar minValue = bounds[axis*2];
        scalar maxValue = bounds[axis*2+1];
        label startingIndex = std::ceil(minValue/labelDistance);
        label finalIndex = std::floor(maxValue/labelDistance);
        label numberOfLabels = finalIndex - startingIndex + 1;

        axesActor_->SetNumberOfLabels(axis, numberOfLabels);

        for (label i = 0; i < numberOfLabels; i++)
        {
            axesActor_->SetLabel(axis, i, static_cast<scalar>(i + startingIndex)*labelDistance);
        }
    }
    else
    {
        axesActor_->SetUseCustomLabels(axis, false);
    }
}

void GridWidget::removeLabels()
{
    if (axesActor_)
    {
        axesActor_->SetUseCustomLabels(X, false);
        axesActor_->SetUseCustomLabels(Y, false);
        axesActor_->SetUseCustomLabels(Z, false);
    }
}



void GridWidget::update(vtkRenderer* renderer)
{
    if (mustUpdateBoundsWithSceneBounds())
    {
        removeFromRenderer(renderer);
        scalar bounds[6];
        boundsUtils::computeAllProcsRenderBoundingBox(renderer, bounds);
        axesActor_->SetGridBounds(bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);
        configureCustomLabels(bounds);
        addToRenderer(renderer);
    }
}

bool GridWidget::mustUpdateBoundsWithSceneBounds() const
{
    return !data_.useCustomBounds;
}



void GridWidget::registerPolyDataToCompositer(vtkRenderer* renderer, Compositer* compositer)
{
    axesActor_->UpdateGeometry(renderer);
    vtkSmartPointer<vtkPropCollection> props = vtkSmartPointer<vtkPropCollection>::New();
    axesActor_->GetActors(props);
    props->InitTraversal();
    for (vtkProp *prop = props->GetNextProp(); prop != nullptr; prop = props->GetNextProp())
    {
        vtkActor* actor = vtkActor::SafeDownCast(prop);
        if (actor)
        {
            compositer->registerWidgetActorPolyData(actor);
        }
    }
}

void GridWidget::redistributePolyDataWithCompositer(Compositer* compositer)
{
    vtkSmartPointer<vtkPropCollection> props = vtkSmartPointer<vtkPropCollection>::New();
    axesActor_->GetActors(props);
    props->InitTraversal();
    for (vtkProp *prop = props->GetNextProp(); prop != nullptr; prop = props->GetNextProp())
    {
        vtkActor* actor = vtkActor::SafeDownCast(prop);
        if (actor)
        {
            compositer->redistributeWidgetActorPolyData(actor);
        }
    }
}


void GridWidget::addToRenderer(vtkRenderer* renderer)
{
    if (axesActor_)
    {
        renderer->AddViewProp(axesActor_);
    }
}

void GridWidget::removeFromRenderer(vtkRenderer* renderer)
{
    if (axesActor_)
    {
        renderer->RemoveViewProp(axesActor_);
    }
}

// ************************************************************************* //
} // End namespace
