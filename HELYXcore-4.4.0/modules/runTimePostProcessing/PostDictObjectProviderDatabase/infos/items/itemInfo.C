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
#include "itemInfo.H"

namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ItemInfo::ItemInfo(const Visualisation* defaultVisualisation) :
    id_(),
    visualisation_(defaultVisualisation)
{
};

ItemInfo::ItemInfo
(
    const dictionary& dict,
    const Visualisation* defaultVisualisation
)
:
    id_(dict),
    visualisation_(dict, defaultVisualisation)
{
}

ItemInfo::ItemInfo
(
    const dictionary& dict,
    ItemType::Value type,
    const Visualisation* defaultVisualisation
)
:
    visualisation_(dict, defaultVisualisation)
{
    id_.readDictTypeless(dict);
    id_.type = ItemType(type);
}

ItemRequirements ItemInfo::getItemRequirements() const
{
    ItemRequirements requirements;
    if (isVisible())
    {
        requirements.getWritableRequiredFields().merge(visualisation_.getRequiredFields());
        bool profile = visualisation_.representation.isProfile();
        bool transparent = !profile && (visualisation_.opacity < 1.0 || visualisation_.representation.isWireframe());
        if ((transparent && needsGhostCellsWhenTransparent()) || (profile && needsGhostCellsWhenProfile()))
        {
            requirements.setNeedsGhostCellsTrue();
        }
    }
    addExtraItemRequirements(requirements);
    return requirements;
}

ItemInfo::~ItemInfo() = default;

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool ItemInfo::isDataEqualTo(const ItemInfo* other) const
{
    if (other == nullptr)
    {
        return false;
    }
    if (other->getId().type != id_.type)
    {
        return false;
    }
    return isSubclassDataEqualTo(other);
}

void ItemInfo::checkObjectType(ItemType::Value expectedType) const
{
    if(id_.type != ItemType(expectedType))
    {
        // This is a code error, not an IO error - this struct should
        // only ever be constructed with the correct type.
        FatalError
            << "Failed to construct object information!"
            << "  Got type \"" << id_.type << "\" from dictionary \""
            << id_.name << "\", but expected "
            << expectedType << exit(FatalError);
    }
}

void ItemInfo::readVisualisationDict(const dictionary& dict)
{
    visualisation_.readFromDict(dict);
}

bool ItemInfo::isVisible() const
{
    return visualisation_.isVisible();
}

bool ItemInfo::isTransparent() const
{
    return visualisation_.isTransparent();
}

const foamField& ItemInfo::getColorField() const
{
    return visualisation_.colourField;
}

}


// ************************************************************************* //
