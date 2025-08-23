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
    (c) 2020 Engys Ltd.

\*---------------------------------------------------------------------------*/

// OpenFOAM includes
#include "colourLookupTable.H"
#include "Utils/Utils.H"
#include "colourMaps/colours.H"

// VTK includes
#include "vtkRenderer.h"
#include "vtkTextProperty.h"
#include "vtkScalarsToColors.h"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimeVis::ColourLookupTable::ColourLookupTable
(
    const ColourLookupTableData& data,
    const ColourMaps& colourMaps
)
:
    data_(data),
    previousDomainMinMax_(-VGREAT, VGREAT)
{
    initializeActorLookupTable();
    initializeLegendLookupTable();
    initializeColorTransferFunction(colourMaps);

    if (!data_.automaticRange)
    {
        buildActorLookupTable(data.fixedRange);
        buildLegendLookupTable(data.fixedRange);
    }
}



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::functionObjects::runTimeVis::ColourLookupTable::updateDomainRange(scalarMinMax domainMinMax)
{
    if (mustUpdateRange(domainMinMax))
    {
        previousDomainMinMax_ = domainMinMax;
        if (!Utils::isRangeValid(domainMinMax) || data_.automaticRange)
        {
            buildActorLookupTable(domainMinMax);
            buildLegendLookupTable(domainMinMax);
        }
        else
        {
            buildActorLookupTable(data_.fixedRange);
            buildLegendLookupTable(data_.fixedRange);
        }
    }
}



void Foam::functionObjects::runTimeVis::ColourLookupTable::initializeColorTransferFunction
(
        const ColourMaps& colourMaps
)
{
    const ColourMap* colourMap = colourMaps.getColourMap(data_.colourMapName);

    colorTransferFunction_ = vtkSmartPointer<vtkDiscretizableColorTransferFunction>::New();
    colorTransferFunction_->DiscretizeOn();
    colorTransferFunction_->SetHSVWrap(colourMap->hsvWrap);
    switch(colourMap->colourSpace.getValue())
    {
    case ColourSpace::Value::RGB:
        colorTransferFunction_->SetColorSpaceToRGB();
        break;
    case ColourSpace::Value::LAB:
        colorTransferFunction_->SetColorSpaceToLab();
        break;
    case ColourSpace::Value::HSV:
        colorTransferFunction_->SetColorSpaceToHSV();
        break;
    case ColourSpace::Value::DIV:
        colorTransferFunction_->SetColorSpaceToDiverging();
        break;
    }

    // Set the ctf control points
    colorTransferFunction_->SetNumberOfValues(data_.resolution);
    for (label i = 0; i < colourMap->colours.size(); ++i)
    {
        label j = i;
        if (data_.inverted)
        {
            j = colourMap->colours.size() - (i+1);
        }

        auto pointX = static_cast<scalar>(colourMap->validPointForIndex(j));

        colorTransferFunction_->AddRGBPoint
        (
            scalar(data_.resolution) * pointX,
            colourMap->colours[i].x(),
            colourMap->colours[i].y(),
            colourMap->colours[i].z()
        );
    }
    colorTransferFunction_->Build();
}

void Foam::functionObjects::runTimeVis::ColourLookupTable::initializeActorLookupTable()
{
    actorLookupTable_ = vtkSmartPointer<vtkLookupTable>::New();
    setSharedLookupTableProperties(actorLookupTable_);
    actorLookupTable_->SetUseAboveRangeColor(0);
    actorLookupTable_->SetUseBelowRangeColor(0);
    actorLookupTable_->SetNumberOfTableValues(data_.resolution);
    actorLookupTable_->SetNanColor(NAN_ACTOR_COLOR[0], NAN_ACTOR_COLOR[1], NAN_ACTOR_COLOR[2], 1);
}

void Foam::functionObjects::runTimeVis::ColourLookupTable::initializeLegendLookupTable()
{
    legendLookupTable_ = vtkSmartPointer<vtkDiscretizableColorTransferFunction>::New();
    setSharedLookupTableProperties(legendLookupTable_);
    legendLookupTable_->SetUseAboveRangeColor(0);
    legendLookupTable_->SetUseBelowRangeColor(0);
    legendLookupTable_->SetNanColor(NAN_BAR_COLOR[0], NAN_BAR_COLOR[1], NAN_BAR_COLOR[2]);
    legendLookupTable_->SetDiscretize(1);
    legendLookupTable_->SetNumberOfValues(data_.resolution);
}

void Foam::functionObjects::runTimeVis::ColourLookupTable::setSharedLookupTableProperties
(
        vtkScalarsToColors* lookupTable
) const
{
    if (data_.fieldName.isComponent())
    {
        lookupTable->SetVectorModeToComponent();
        lookupTable->SetVectorComponent(data_.fieldName.getComponentIndex());
    }
    else  // Either a scalar field or a vector magnitude field
    {
        lookupTable->SetVectorModeToMagnitude();
    }
}

Foam::scalarMinMax Foam::functionObjects::runTimeVis::ColourLookupTable::adjustRange(scalarMinMax range) {
    static const scalar ABSOLUTE_ADJUSTMENT = std::pow(Foam::VSMALL, static_cast<scalar>(2.0) / static_cast<scalar>(3.0));
    static const scalar RELATIVE_ADJUSTMENT = 1 + std::pow(Foam::SMALL, static_cast<scalar>(2.0) / static_cast<scalar>(3.0));
    scalar start, end, minEnd;
    start = range.min();
    end = range.max();
    minEnd = std::nextafter(start + ABSOLUTE_ADJUSTMENT, start + static_cast<scalar>(1.0)) * RELATIVE_ADJUSTMENT;
    end = std::max(minEnd, end);
    return {start, end};
}

void Foam::functionObjects::runTimeVis::ColourLookupTable::buildActorLookupTable(scalarMinMax range)
{
    if (Utils::isRangeValid(range))
    {
        actorLookupTable_->SetTableRange(range.min(), range.max());

        double c[3];
        for (label i = 0; i < data_.resolution; i++)
        {
            colorTransferFunction_->GetColor(i, c);
            actorLookupTable_->SetTableValue(i, c[0], c[1], c[2], 1.0);
        }
    }
    else
    {
        actorLookupTable_->SetTableRange(0, Foam::VSMALL);
    }
    actorLookupTable_->Build();
}

void Foam::functionObjects::runTimeVis::ColourLookupTable::buildLegendLookupTable(scalarMinMax range)
{
    legendLookupTable_->RemoveAllPoints();

    if (Utils::isRangeValid(range))
    {
        scalarMinMax adjustedRange = adjustRange(range);
        setColorTableProperties(adjustedRange);
    }
    else
    {
        setNanTableProperties();
    }
    legendLookupTable_->Build();
}

void Foam::functionObjects::runTimeVis::ColourLookupTable::setColorTableProperties(Foam::scalarMinMax &adjustedRange)
{
    if (data_.resolution <= 1)
    {
        if (data_.resolution == 1)
        {
            double c[3];
            colorTransferFunction_->GetColor(0, c);
            legendLookupTable_->AddRGBPoint(adjustedRange.min(), c[0], c[1], c[2]);
            legendLookupTable_->AddRGBPoint(adjustedRange.max(), c[0], c[1], c[2]);
        }
        else
        {
            legendLookupTable_->AddRGBPoint(adjustedRange.min(), 0, 0, 0);
            legendLookupTable_->AddRGBPoint(adjustedRange.max(), 0, 0, 0);
        }
    } else
    {
        std::vector<scalar> indexes = Utils::linspace(adjustedRange.min(), adjustedRange.max(), data_.resolution);
        for (label i = 0; i < data_.resolution; i++)
        {
            double c[3];
            colorTransferFunction_->GetColor(i, c);
            legendLookupTable_->AddRGBPoint(indexes[i], c[0], c[1], c[2]);
        }
    }
}

void Foam::functionObjects::runTimeVis::ColourLookupTable::setNanTableProperties()
{
    legendLookupTable_->SetNanColor(NAN_BAR_COLOR[0], NAN_BAR_COLOR[1], NAN_BAR_COLOR[2]);
    legendLookupTable_->SetNumberOfValues(2);
    legendLookupTable_->AddRGBPoint(0, NAN_BAR_COLOR[0], NAN_BAR_COLOR[1], NAN_BAR_COLOR[2]);
    legendLookupTable_->AddRGBPoint(1e-100, NAN_BAR_COLOR[0], NAN_BAR_COLOR[1], NAN_BAR_COLOR[2]);
}

bool Foam::functionObjects::runTimeVis::ColourLookupTable::mustUpdateRange(scalarMinMax domainMinMax)
{
    if (Utils::isRangeValid(domainMinMax) != Utils::isRangeValid(previousDomainMinMax_))
    {
        return true;
    }
    else if (data_.automaticRange)
    {
        return (domainMinMax.min() != previousDomainMinMax_.min() ||
                domainMinMax.max() != previousDomainMinMax_.max());
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
