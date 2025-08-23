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
    (c) 2020-2021 Engys Ltd.

\*---------------------------------------------------------------------------*/

// OpenFOAM includes
#include "colorLegends.H"

#include "colourLookupTable/colourLookupTableProvider.H"
#include "storage/foamMeshes.H"

// VTK includes
#include "vtkRenderer.h"

namespace Foam::functionObjects::runTimeVis
{

std::unique_ptr<ColorLegends> ColorLegends::createWidgetIfVisible(const ColourLegendsData& data)
{
    if (data.isAnyColorLegendVisible())
    {
        return std::unique_ptr<ColorLegends>(new ColorLegends(data));
    }
    return {nullptr};
}

ColorLegends::ColorLegends(const ColourLegendsData& data)
{
    for (const ColourLegendData& colourLegendData : data.colourLegends)
    {
        ColorLegend colorLegend(colourLegendData);
        if (colorLegend.isVisible())
        {
            colourLegendScalarBars_.insert(colourLegendData.name, colorLegend);
        }
    }
}

void ColorLegends::update(
        ColourLookupTableProvider& colourLutProvider,
        const ColourMaps& colourMaps,
        const FoamMeshes& meshes,
        const ExternalFields& externalFields
)
{
    for (ColorLegend& s : colourLegendScalarBars_)
    {
        scalarMinMax domainMinMax = meshes.getDomainRangeForField(s.getName());
        domainMinMax += externalFields.getDomainRangeForField(s.getName());
        const ColourLookupTable* lut = colourLutProvider.updateAndReturnColorLookupTable(s.getName(), colourMaps, domainMinMax);
        s.updateLookupTable(lut->getLegendLookupTable(), domainMinMax);
    }
}


void ColorLegends::addToRenderer(vtkRenderer* renderer)
{
    for (ColorLegend& s : colourLegendScalarBars_)
    {
        s.addToRenderer(renderer);
    }
}

void ColorLegends::removeFromRenderer(vtkRenderer* renderer) {
    for (ColorLegend& s : colourLegendScalarBars_)
    {
        s.removeFromRenderer(renderer);
    }
}


} // End namespace Foam

// ************************************************************************* //
