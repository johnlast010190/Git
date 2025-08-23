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
    (c) 2022-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

// OpenFOAM includes
#include "sceneInfo.H"
#include "containers/Lists/DynamicList/DynamicList.H"

namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SceneInfo::SceneInfo
(
    const dictionary& sceneDict,
    ItemInfoFactory& factory,
    const ReferenceFrames& rfs
)
:
    cameras(rfs, sceneDict),
    backgroundColours(GET_OPTIONAL_DICTIONARY(sceneDict,sceneKeys::BACKGROUND_COLORS_DICT_KEY)),
    widgetsInfo(sceneDict),
    colourLookupTablesInfo(sceneDict)
{
    name = Utils::getSubdictName(sceneDict);

    auto itemsDictList = sceneDict.lookup<List<dictionary>>(sceneKeys::ITEMS_DICT_KEY);
    addItems(itemsDictList, factory);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void SceneInfo::addItems(const List<dictionary>& itemsDictList, ItemInfoFactory& factory)
{
    itemInfos.clear();

    for(const dictionary& itemDict : itemsDictList)
    {
        const ItemInfo* sceneItemInfo = factory.createAndReturnItemInfoForSceneItem(itemDict, name);
        itemInfos.insert(sceneItemInfo->getId(), sceneItemInfo);
    }
}

ItemRequirements SceneInfo::getItemRequirements() const
{
    ItemRequirements fields;
    for (const ItemInfo* itemInfo : itemInfos)
    {
        fields.mergeFromDownstream(itemInfo->getItemRequirements());
    }
    fields.mergeFromDownstream(widgetsInfo.colourLegendsData.getItemRequirements());
    return fields;
}



// ************************************************************************* //
}
