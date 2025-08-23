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

Class
    Foam::functionObjects::runTimeVis::SceneInfo

Description
    Information about a scene
    Contains the data from the dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "baseClasses/itemRequirements.H"

#include "infos/widgetsInfo.H"
#include "infos/colourLookupTablesInfo.H"

#include "infos/items/itemInfo.H"
#include "infos/items/itemInfoFactory.H"

#include "dataStructs/scene/backgroundColourData.H"
#include "dataStructs/scene/camerasData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct SceneInfo
{
    std::string name;
    CamerasData cameras;
    BackgroundColourData backgroundColours;

    WidgetsInfo widgetsInfo;
    ColourLookupTablesInfo colourLookupTablesInfo;

    HashTable<const ItemInfo*, Id, IdHasher> itemInfos;

    // Construct from dictionary
    SceneInfo(const dictionary& sceneDict, ItemInfoFactory& factory, const ReferenceFrames& rfs);

    [[nodiscard]] ItemRequirements getItemRequirements() const;

private:

    void addItems(const List<dictionary>& itemsDictList, ItemInfoFactory& factory);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
