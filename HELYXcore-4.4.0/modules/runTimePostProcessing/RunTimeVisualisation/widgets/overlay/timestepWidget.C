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

#include "timestepWidget.H"

#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkRenderer.h"

namespace Foam::functionObjects::runTimeVis
{

std::unique_ptr<TimestepWidget> TimestepWidget::createWidgetIfVisible(const TimestepWidgetData& data)
{
    if (data.visible)
    {
        return std::unique_ptr<TimestepWidget>(new TimestepWidget(data));
    }
    else
    {
        return nullptr;
    }
}

TimestepWidget::TimestepWidget(const TimestepWidgetData& data)
:
        data_(data)
{
    timestepWidget_ = vtkSmartPointer<vtkTextActor>::New();

    timestepWidget_->SetTextScaleModeToNone();

    vtkTextProperty* textProperty = timestepWidget_->GetTextProperty();
    textProperty->SetColor(data.font.colour[0], data.font.colour[1], data.font.colour[2]);
    
    if (data.font.type.IsCourier())
    {
        textProperty->SetFontFamilyToCourier();
    }
    else if (data.font.type.IsTimes())
    {
        textProperty->SetFontFamilyToTimes();
    }
    else
    {
        textProperty->SetFontFamilyToArial();
    }

    textProperty->SetFontSize(data.font.size);
    textProperty->SetJustificationToRight();
    textProperty->SetVerticalJustificationToTop();
    textProperty->SetOpacity(data.font.opacity);
    textProperty->SetBold(data.font.bold ? 1 : 0);
    textProperty->SetItalic(data.font.italics ? 1 : 0);
    textProperty->SetShadow(data.font.shadow ? 1 : 0);

    timestepWidget_->GetPositionCoordinate()->SetCoordinateSystemToNormalizedDisplay();
    timestepWidget_->GetPosition2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
    timestepWidget_->GetPositionCoordinate()->SetValue
    (
        data.coordinates.x(),
        data.coordinates.y()
    );

    timestepWidget_->SetInput("");
}



void TimestepWidget::update(scalar currentTimeValue)
{
    int size_s = std::snprintf( nullptr, 0, data_.format.c_str(), currentTimeValue) + 1;
    if (size_s <= 0) {
        Warning << "Invalid timestep widget format. Using default.";
        timestepWidget_->SetInput(std::to_string(currentTimeValue).c_str());
    } else {
        std::vector<char> buffer(size_s + 1);
        std::snprintf(buffer.data(), size_s, data_.format.c_str(), currentTimeValue);
        buffer[size_s] = 0;
        timestepWidget_->SetInput(buffer.data());
    }
}

void TimestepWidget::addToRenderer(vtkRenderer* renderer)
{
    renderer->AddActor2D(timestepWidget_);
}

void TimestepWidget::removeFromRenderer(vtkRenderer* renderer)
{
    renderer->RemoveActor2D(timestepWidget_);
}


// ************************************************************************* //
} // End namespace
