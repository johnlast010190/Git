/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
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
    (c) 2020-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

// OpenFOAM includes
#include <utility>

#include "items.H"

#include "postDictObjectProviderDatabase.H"
#include "Utils/vtkDataHandlingTools.H"

// VTK includes
#include "vtkActor.h"
#include "engysAppendFilter.h"
#include "vtkPointData.h"

namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Items::Items
    (
        const PostDictObjectProviderDatabase &database,
        const HashTable<const ItemInfo *, Id, IdHasher> &itemInfos,
        std::string sceneName
    )
    :
    database_(database),
    sceneName_(std::move(sceneName))
{
    items_.clear();
    int index = 0;
    for (const ItemInfo *itemInfo: itemInfos)
    {
        if (itemInfo->isVisible())
        {
            items_.emplace_back(itemInfo, index);
        }
        if (itemInfo->getColorField().isIndexedColour())
        {
            index++;
        }
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Items::updateActors
    (
        label currentTimeLabel,
        scalar currentTime,
        ColourLookupTableProvider &provider,
        const ColourMaps &colourMaps
    )
{
    // update and colour each actor
    for (Item &item: items_)
    {
        rtppActor *actor = item.getConfiguredActor();

        if (actor->isOutline())
        {
            vtkSmartPointer<vtkPolyData> outline = database_.getOutlineForSceneItem(
                item.getInfo(),
                sceneName_,
                currentTimeLabel,
                currentTime
            );
            item.setActorOutline(outline);
        }
        else
        {
            vtkSmartPointer<vtkDataSet> data = database_.getDataSetForSceneItem(
                item.getInfo(),
                sceneName_,
                currentTimeLabel,
                currentTime
            );
            item.setActorDataSet(data);
        }

        scalarMinMax range = database_.getFoamMeshes().getDomainRangeForField(item.getColorField());
        range += database_.getExternalFields().getDomainRangeForField(item.getColorField());
        if (!item.getColorField().isSolidColour())
        {
            const ColourLookupTable *lut = provider.updateAndReturnColorLookupTable(
                item.getColorField(),
                colourMaps,
                range
            );
            actor->updateLookupTable(lut->getActorLookupTable());
        }
    }
}

vtkSmartPointer<vtkDataSet> Items::createDataSetWithAllData
    (
        label currentTimeLabel,
        scalar currentTime
    )
{
    vtkNew<engysAppendFilter> appendFilter;
    for (Item &item: items_)
    {
        vtkSmartPointer<vtkDataSet> data = database_.getDataSetForSceneItem(
            item.getInfo(),
            sceneName_,
            currentTimeLabel,
            currentTime
        );
        if (data->GetNumberOfPoints() <= 0 && data->GetNumberOfCells() <= 0)
        {
            continue;
        }
        appendFilter->AddInputData(vtk::Tools::removeGhostCells(data));
    }

    appendFilter->Update();
    vtkSmartPointer<vtkDataSet> output = appendFilter->GetOutput();
    output->GetPointData()->RemoveArray("engysGhostMeshPointsArray");
    return output;
}

void Items::registerDataToCompositer(Compositer *compositer)
{
    for (Item &item: items_)
    {
        item.registerActorDataSetToCompositer(compositer);
    }
}

void Items::redistributeDataWithCompositerIfNecessary(Compositer *compositer)
{
    for (Item &item: items_)
    {
        item.redistributeActorDataSetIfNecessary(compositer);
    }
}

void Items::addToRenderer(renderManager &renderer)
{
    for (Item &item: items_)
    {
        rtppActor *actor = item.getConfiguredActor();
        if (actor)
        {
            renderer.addActor(*actor);
        }
    }
}

bool Items::isAnyVisible() const
{
    return std::any_of(items_.begin(), items_.end(), [](const Item &item) { return item.isVisible(); });
}


} // End namespace Foam

// ************************************************************************* //
