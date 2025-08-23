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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

// OpenFOAM includes
#include "widgetsInOverlay.H"
#include "infos/widgetsInfo.H"
#include "colourLookupTable/colourLookupTableProvider.H"
#include "widgets/overlay/axesWidget.H"
#include "widgets/overlay/logoWidget.H"
#include "widgets/overlay/timestepWidget.H"
#include "widgets/overlay/colorLegends.H"
#include "widgets/overlay/vectorWidget.H"
#include "Utils/ParallelUtils.H"

#include "vtkRenderer.h"

#define callOverlayWidgetMethodIfVisible(widget, method) \
    if ((widget).get()) \
    { \
        (widget)->method; \
    } \

#define createWidgetIfVisibleAndMaster(widgetClass, ...) \
    ParallelUtils::isMaster() ? widgetClass::createWidgetIfVisible(__VA_ARGS__) : nullptr

namespace Foam::functionObjects::runTimeVis
{

WidgetsInOverlay::WidgetsInOverlay(const WidgetsInfo& widgetsInfo)
:
        renderer_(nullptr),
        axesWidget_(createWidgetIfVisibleAndMaster(AxesWidget, widgetsInfo.axisWidgetData)),
        colorLegends_(createWidgetIfVisibleAndMaster(ColorLegends, widgetsInfo.colourLegendsData)),
        logoWidget_(createWidgetIfVisibleAndMaster(LogoWidget, widgetsInfo.logoWidgetData)),
        timestepWidget_(createWidgetIfVisibleAndMaster(TimestepWidget, widgetsInfo.timestepWidgetData)),
        vectorWidget_(createWidgetIfVisibleAndMaster(VectorWidget, widgetsInfo.vectorWidgetData))
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void WidgetsInOverlay::addVisibleWidgetsToRenderer()
{
    callOverlayWidgetMethodIfVisible(axesWidget_, addToRenderer(renderer_))
    callOverlayWidgetMethodIfVisible(colorLegends_, addToRenderer(renderer_))
    callOverlayWidgetMethodIfVisible(logoWidget_, addToRenderer(renderer_))
    callOverlayWidgetMethodIfVisible(timestepWidget_, addToRenderer(renderer_))
    callOverlayWidgetMethodIfVisible(vectorWidget_, addToRenderer(renderer_))
}

void WidgetsInOverlay::updateVisibleWidgets(
        scalar currentTimeValue,
        ColourLookupTableProvider& colourLutProvider,
        const ColourMaps& colourMaps,
        const FoamMeshes& meshes,
        const ExternalFields& externalFields,
        const Time& runTime
)
{
    callOverlayWidgetMethodIfVisible(timestepWidget_, update(currentTimeValue))
    callOverlayWidgetMethodIfVisible(colorLegends_, update(colourLutProvider, colourMaps, meshes, externalFields))
    callOverlayWidgetMethodIfVisible(vectorWidget_, update(runTime))
}


} // End namespace

// ************************************************************************* //
