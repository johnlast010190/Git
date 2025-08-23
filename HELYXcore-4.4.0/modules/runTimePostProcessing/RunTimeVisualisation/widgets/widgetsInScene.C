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
#include "widgetsInScene.H"
#include "infos/widgetsInfo.H"
#include "storage/referenceFrames.H"

#define callSceneWidgetMethodIfVisible(widget, method) \
    if ((widget).isVisible()) \
    { \
        (widget).method; \
    } \


namespace Foam::functionObjects::runTimeVis
{

WidgetsInScene::WidgetsInScene(Compositer* compositer, const WidgetsInfo& widgetsInfo, const ReferenceFrames& referenceFrames)
:
        baseRenderer_(nullptr),
        compositer_(compositer),
        gridWidget_(widgetsInfo.gridWidgetData)
{
    for (const word& referenceFrameName : widgetsInfo.referenceFramesWidgetData.referenceFrames)
    {
        referenceFrameWidgets_.emplace_back(referenceFrames.getReferenceFrame(referenceFrameName));
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void WidgetsInScene::addVisibleWidgetsToRenderer()
{
    callSceneWidgetMethodIfVisible(gridWidget_, addToRenderer(baseRenderer_))
    for (const ReferenceFrameWidget& referenceFrameWidget : referenceFrameWidgets_)
    {
        referenceFrameWidget.addToRenderer(baseRenderer_);
    }
}

void WidgetsInScene::updateVisibleWidgets()
{
    callSceneWidgetMethodIfVisible(gridWidget_, update(baseRenderer_))
    for (ReferenceFrameWidget& referenceFrameWidget : referenceFrameWidgets_)
    {
        referenceFrameWidget.update();
    }
}

void WidgetsInScene::registerDataToCompositer()
{
    callSceneWidgetMethodIfVisible(gridWidget_, registerPolyDataToCompositer(baseRenderer_, compositer_))
    for (const ReferenceFrameWidget& referenceFrameWidget : referenceFrameWidgets_)
    {
        referenceFrameWidget.registerPolyDataToCompositer(compositer_);
    }
}

void WidgetsInScene::redistributeDataWithCompositer()
{
    callSceneWidgetMethodIfVisible(gridWidget_, redistributePolyDataWithCompositer(compositer_))
    for (const ReferenceFrameWidget& referenceFrameWidget : referenceFrameWidgets_)
    {
        referenceFrameWidget.redistributePolyDataWithCompositer(compositer_);
    }
}


} // End namespace

// ************************************************************************* //
