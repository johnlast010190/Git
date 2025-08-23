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
    (c) 2019 OpenCFD Ltd.
    (c) 2020-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

// OpenFOAM includes
#include "colorLegend.H"
#include "primitives/Scalar/scalar/scalar.H"
#include "Utils/Utils.H"

// VTK includes
#include "vtkDiscretizableColorTransferFunction.h"
#include "vtkRenderer.h"
#include "vtkSmartPointer.h"
#include "vtkTextProperty.h"
#include "vtkContext2DScalarBarRepresentation.h"
#include "vtkScalarBarWidget.h"

#define AUTO_RANGE_LABEL_E_FORMAT "%6.2e"
#define AUTO_RANGE_LABEL_F_FORMAT "%6.2f"
#define POSITION_MARGIN 0.03

#if defined(HELYX_SP)
#define ABS fabsf
#else
#define ABS fabs
#endif

namespace Foam::functionObjects::runTimeVis
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ColorLegend::ColorLegend
    (
        const ColourLegendData &legendData
    )
    :
    legendData_(legendData)
{
    initialiseActor();
}


void ColorLegend::updateLookupTable
    (
        vtkDiscretizableColorTransferFunction *lut,
        const scalarMinMax &domainRange
    )
{
    static const scalar RELATIVE_TOL = std::pow(Foam::SMALL, static_cast<scalar>(2.0)/static_cast<scalar>(3.0));
    static const scalar ABSOLUTE_TOL = std::pow(Foam::ROOTVSMALL, static_cast<scalar>(2.0)/static_cast<scalar>(3.0));

    double *tempRange = lut->GetRange();
    Foam::scalarMinMax legendRange(static_cast<scalar>(tempRange[0]), static_cast<scalar>(tempRange[1]));
    actor_->SetLookupTable(lut);

    bool autoLabelFormat = legendData_.autoLabels;
    scalar maxAbs = std::max(ABS(legendRange.max()), ABS(legendRange.min()));
    scalar diffAbs = ABS(legendRange.max() - legendRange.min());
    if (autoLabelFormat && (diffAbs <= maxAbs * RELATIVE_TOL || diffAbs <= ABSOLUTE_TOL))
    {
        // if the legend is set to auto label and the difference between the ranges is too small, then force this format
        actor_->SetAutomaticLabelFormat(0);
        actor_->SetLabelFormat(AUTO_RANGE_LABEL_F_FORMAT);
        actor_->SetRangeLabelFormat(AUTO_RANGE_LABEL_F_FORMAT);
    }
    else
    {
        actor_->SetAutomaticLabelFormat(autoLabelFormat ? 1 : 0);
        actor_->SetLabelFormat(legendData_.labelFormat.c_str());
        if (autoLabelFormat)
        {
            actor_->SetRangeLabelFormat(decideAutoRangesFormat(legendRange.min(), legendRange.max()).c_str());
        }
        else
        {
            actor_->SetRangeLabelFormat(legendData_.labelFormat.c_str());
        }
    }

    if (Utils::isRangeValid(domainRange))
    {
        actor_->SetUseCustomLabels(!autoLabelFormat);
        label numberOfLabels = legendData_.numberOfLabels;
        if (numberOfLabels >= 1)
        {
            std::vector<Foam::scalar> values = Utils::linspace(legendRange.min(), legendRange.max(), numberOfLabels);
            for (label i = 0; i < numberOfLabels; i++)
            {
                actor_->SetCustomLabel(i, values.at(i));
            }
        }
    }
    else
    {
        // if the legend is set to auto label and range is invalid
        actor_->SetUseCustomLabels(true);
        actor_->SetAddRangeLabels(0);
        actor_->SetNumberOfCustomLabels(0);
    }

    updatePositioning(representation_);
}


void ColorLegend::initialiseActor()
{
    actor_ = vtkSmartPointer<vtkContext2DScalarBarActor>::New();
    actor_->GetPositionCoordinate()->SetCoordinateSystemToNormalizedDisplay();
    actor_->SetVerticalTitleSeparation(10);
    actor_->VisibilityOff();

    actor_->GetLabelTextProperty()->SetColor
        (
            legendData_.font.colour[0],
            legendData_.font.colour[1],
            legendData_.font.colour[2]
        );
    actor_->GetLabelTextProperty()->SetFontSize(legendData_.font.size);
    actor_->GetLabelTextProperty()->SetOpacity(legendData_.font.opacity);
    actor_->GetTitleTextProperty()->SetColor
        (
            legendData_.font.colour[0],
            legendData_.font.colour[1],
            legendData_.font.colour[2]
        );
    actor_->GetTitleTextProperty()->SetFontSize(legendData_.font.size);
    actor_->GetTitleTextProperty()->SetOpacity(legendData_.font.opacity);

    actor_->GetLabelTextProperty()->SetBold(legendData_.font.bold ? 1 : 0);
    actor_->GetLabelTextProperty()->SetItalic(legendData_.font.italics ? 1 : 0);
    actor_->GetLabelTextProperty()->SetShadow(legendData_.font.shadow ? 1 : 0);
    actor_->GetTitleTextProperty()->SetBold(legendData_.font.bold ? 1 : 0);
    actor_->GetTitleTextProperty()->SetItalic(legendData_.font.italics ? 1 : 0);
    actor_->GetTitleTextProperty()->SetShadow(legendData_.font.shadow ? 1 : 0);

    actor_->SetTitle(stringOps::center(legendData_.title, 10).c_str());

    if (legendData_.font.type.IsCourier())
    {
        actor_->GetLabelTextProperty()->SetFontFamilyToCourier();
        actor_->GetTitleTextProperty()->SetFontFamilyToCourier();
    }
    else if (legendData_.font.type.IsTimes())
    {
        actor_->GetLabelTextProperty()->SetFontFamilyToTimes();
        actor_->GetTitleTextProperty()->SetFontFamilyToTimes();
    }
    else
    {
        actor_->GetLabelTextProperty()->SetFontFamilyToArial();
        actor_->GetTitleTextProperty()->SetFontFamilyToArial();
    }

    bool autoLabelFormat = legendData_.autoLabels;

    actor_->SetUseCustomLabels(!autoLabelFormat);
    label numberOfLabels = legendData_.numberOfLabels;
    if (numberOfLabels >= 1)
    {
        if (numberOfLabels >= 2)
        {
            actor_->SetAddRangeLabels(1);
        }
        else
        {
            actor_->SetAddRangeLabels(0);
        }
        actor_->SetNumberOfCustomLabels(numberOfLabels);
    }
    else
    {
        actor_->SetAddRangeLabels(0);
        actor_->SetNumberOfCustomLabels(0);
    }

    actor_->SetScalarBarLength(legendData_.length);
    actor_->SetScalarBarThickness(legendData_.thickness);
    actor_->ForceHorizontalTitleOn();
    actor_->SetDrawTickMarks(legendData_.showTicks);
    actor_->SetTitleJustification(VTK_TEXT_CENTERED);

    representation_ = vtkSmartPointer<vtkContext2DScalarBarRepresentation>::New();
    representation_->SetScalarBarActor(actor_);
    representation_->SetShowBorderToOff();
    representation_->SetBorderMargin(POSITION_MARGIN);

    if (legendData_.visible)
    {
        actor_->VisibilityOn();
    }
}

void ColorLegend::addToRenderer(vtkRenderer* renderer) const
{
    renderer->AddActor2D(representation_);
}

void ColorLegend::removeFromRenderer(vtkRenderer* renderer) const
{
    renderer->RemoveActor2D(representation_);
}

bool ColorLegend::isVertical() const
{
    if (legendData_.location.IsCustom())
    {
        return legendData_.vertical;
    }
    else if (legendData_.location.IsCenterX())
    {
        return false;
    }
    else
    {
        return true;
    }
}

Foam::string ColorLegend::decideAutoRangesFormat(scalar min, scalar max)
{
    if (shouldUseScientificNotation(min) || shouldUseScientificNotation(max))
    {
        return AUTO_RANGE_LABEL_E_FORMAT;
    }
    else
    {
        return AUTO_RANGE_LABEL_F_FORMAT;
    }
}

bool ColorLegend::shouldUseScientificNotation(scalar range)
{
    return (range > 1e4) ||
           (range < -1e4) ||
           (range != 0 && range < 1e-1 && range > -1e-1);
}

void ColorLegend::updatePositioning(vtkContext2DScalarBarRepresentation* representation) const
{
    representation->SetOrientation(isVertical() ? 1 : 0);

    vtkContext2DScalarBarRepresentation::PositionType vtkPosition;
    if (legendData_.location.IsCustom())
    {
        vtkPosition = vtkContext2DScalarBarRepresentation::CUSTOM;
        representation->SetPosition(legendData_.coordinates.first(), legendData_.coordinates.second());
    }
    else
    {
        if (legendData_.location.IsLeft())
        {
            if (legendData_.location.IsTop())
            {
                vtkPosition = vtkContext2DScalarBarRepresentation::LEFT_TOP;
            }
            else if (legendData_.location.IsBottom())
            {
                vtkPosition = vtkContext2DScalarBarRepresentation::LEFT_BOTTOM;
            }
            else
            {
                vtkPosition = vtkContext2DScalarBarRepresentation::LEFT;
            }
        }
        else if (legendData_.location.IsRight())
        {
            if (legendData_.location.IsTop())
            {
                vtkPosition = vtkContext2DScalarBarRepresentation::RIGHT_TOP;
            }
            else if (legendData_.location.IsBottom())
            {
                vtkPosition = vtkContext2DScalarBarRepresentation::RIGHT_BOTTOM;
            }
            else
            {
                vtkPosition = vtkContext2DScalarBarRepresentation::RIGHT;
            }
        }
        else
        {
            if (legendData_.location.IsTop())
            {
                vtkPosition = vtkContext2DScalarBarRepresentation::TOP;
            }
            else
            {
                vtkPosition = vtkContext2DScalarBarRepresentation::BOTTOM;
            }
        }
    }
    representation->SetPositioning(vtkPosition);
}

// ************************************************************************* //
}