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

#include "referenceFrameWidget.H"

#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "engysCoordAxesActor.h"

namespace Foam::functionObjects::runTimeVis
{
ReferenceFrameWidget::ReferenceFrameWidget(const ReferenceFrame &referenceFrame)
    :
    referenceFrame_(referenceFrame)
{
    actor_ = vtkSmartPointer<engysCoordAxesActor>::New();

    actor_->UseFixedScreenSizeOn();
    actor_->SetEditTypeToNone();
    actor_->SetScreenSize(0.1);
    actor_->VisibilityOn();
}

void ReferenceFrameWidget::addToRenderer(vtkRenderer* renderer) const
{
    renderer->AddActor(actor_);
}

void ReferenceFrameWidget::removeFromRenderer(vtkRenderer* renderer) const
{
    renderer->RemoveActor(actor_);
}

void ReferenceFrameWidget::update()
{
    vector secondVector;
    switch(referenceFrame_.getType().getValue())
    {
        case CoordinateSystemType::CARTESIAN:
            actor_->SetCoordinateSystemType(engysCoordAxesActor::CoordinateSystemType::CARTESIAN);
            secondVector = referenceFrame_.e2();
            break;
        case CoordinateSystemType::CYLINDRICAL:
            secondVector = referenceFrame_.e3();
            actor_->SetCoordinateSystemType(engysCoordAxesActor::CoordinateSystemType::CYLINDRICAL);
            break;
        default:
            actor_->VisibilityOff();
            Warning << "Unknown reference frame type: " << referenceFrame_.getType().getName() << endl;
            return;
    }
    point origin = referenceFrame_.getOrigin();
    actor_->SetOrigin(origin.x(), origin.y(), origin.z());

    vector firstVector = referenceFrame_.e1();
    actor_->SetPrimaryAxis(firstVector.x(), firstVector.y(), firstVector.z());

    actor_->SetSecondaryAxis(secondVector.x(), secondVector.y(), secondVector.z());
}

void ReferenceFrameWidget::registerPolyDataToCompositer(Compositer* compositer) const
{
    vtkSmartPointer<vtkPropCollection> props = vtkSmartPointer<vtkPropCollection>::New();
    actor_->GetActors(props);
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

void ReferenceFrameWidget::redistributePolyDataWithCompositer(Compositer* compositer) const
{
    vtkSmartPointer<vtkPropCollection> props = vtkSmartPointer<vtkPropCollection>::New();
    actor_->GetActors(props);
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

// ************************************************************************* //
} // End namespace Foam
