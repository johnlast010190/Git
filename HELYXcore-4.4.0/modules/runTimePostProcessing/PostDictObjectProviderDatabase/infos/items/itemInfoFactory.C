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
    (c) 2022-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include <memory>

#include "itemInfoFactory.H"

#include "infos/items/invalid/missingGeometryInfo.H"
#include "infos/items/invalid/missingItemInfo.H"
#include "infos/items/objects/sliceInfo.H"
#include "infos/items/objects/clipInfo.H"
#include "infos/items/objects/transformInfo.H"
#include "infos/items/objects/isosurfaceInfo.H"
#include "infos/items/objects/streamlinesInfo.H"
#include "infos/items/objects/glyphInfo.H"
#include "infos/items/objects/thresholdInfo.H"
#include "infos/items/objects/fieldSamplingInfo.H"
#include "infos/items/objects/turboPost/turboSliceStreamwiseInfo.H"
#include "infos/items/objects/turboPost/turboSliceSpanwiseInfo.H"
#include "infos/items/objects/turboPost/meridionalGridInfo.H"
#include "infos/items/objects/turboPost/bladeLoadingLineInfo.H"
#include "infos/items/objects/turboPost/hubToShroudLineInfo.H"
#include "infos/items/objects/turboPost/inletToOutletLineInfo.H"
#include "infos/items/objects/importDataInfo.H"
#include "infos/items/testObjects/surfaceSplitterInfo.H"
#include "infos/items/mesh/cellZoneInfo.H"
#include "infos/items/mesh/faceZoneInfo.H"
#include "infos/items/mesh/internalBoundaryInfo.H"
#include "infos/items/mesh/patchInfo.H"
#include "infos/items/mesh/processorBoundaryInfo.H"
#include "infos/items/mesh/volMeshInfo.H"
#include "infos/items/objects/line3DInfo.H"
#include "infos/items/group/groupInfo.H"
#include "infos/items/geometry/boundarySurfaceInfo.H"
#include "infos/items/geometry/surfaceInfos.H"
#include "infos/items/geometry/featureLineInfo.H"
#include "infos/items/fileSource/pvdFileInfo.H"
#include "infos/items/fileSource/fieldSamplingFunctionObjectInfo.H"
#include "infos/items/testObjects/tool3Dto2DInfo.H"
#include "postDict/postDictKeys.H"

#include "types/surfaceType.H"

#include "primitives/strings/fileName/fileName.H"

namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * *  Public Member Functions  * * * * * * * * * * * //

const ItemInfo* ItemInfoFactory::createAndReturnItemInfoForSceneItem(const dictionary &dict, const word& sceneName)
{
    return createStoreAndReturnItemInfoForItem(
        getSceneObjectItemDict(dict),
        dict.lookup<string>(idKeys::NAME_KEY).c_str(),
        dict,
        dict.lookupOrDefault<word>(idKeys::REGION_KEY, regionKeys::DEFAULT_REGION_KEY).c_str(),
        dict.lookup<ItemType>(idKeys::TYPE_KEY).getValue(),
        sceneName.c_str());
}

const dictionary& ItemInfoFactory::getSceneObjectItemDict(const dictionary &dict)
{
    static const dictionary empty;
    auto itemType = dict.lookup<ItemType>(idKeys::TYPE_KEY);
    if (itemType.getValue() == ItemType::SURFACE)
    {
        if (dict.isDict(surfaceKeys::SURFACE_TYPE_SUBDICT_KEY))
        {
            const dictionary& subdict = dict.subDict(surfaceKeys::SURFACE_TYPE_SUBDICT_KEY);
            if (subdict.toc().size() == 1)
            {
                return subdict.subDict(subdict.toc()[0]);
            }
            FatalError << "When defining a new surface inside a scene, it must contain a "
                       << surfaceKeys::SURFACE_TYPE_SUBDICT_KEY << " sub dictionary that contains "
                       << "a sub dictionary named with the file name." << endl
                       << "Error for item" << dict.name() << endl
                       << abort(FatalError);
        }

        return empty;
    }
    else
    {
        return dict;
    }
}

void ItemInfoFactory::createItemInfoForObject
    (
        const dictionary &dict,
        const dictionary &objectVisDict
    )
{
    createStoreAndReturnItemInfoForItem(dict, nullptr, objectVisDict, nullptr, -1, nullptr);
}

void ItemInfoFactory::createItemInfoForPatch
    (
        const dictionary &dict,
        const string &region,
        const string &patchName
    )
{
    createStoreAndReturnItemInfoForItem(dict, patchName.c_str(), dict, region.c_str(), ItemType::PATCH, nullptr);
}

void ItemInfoFactory::createItemInfoForCellZone
    (
        const dictionary &dict,
        const string &region,
        const string &patchName
    )
{
    createStoreAndReturnItemInfoForItem(dict, patchName.c_str(), dict, region.c_str(), ItemType::CELLZONE, nullptr);
}

void ItemInfoFactory::createItemInfoForFaceZone
    (
        const dictionary &dict,
        const string &region,
        const string &patchName
    )
{
    createStoreAndReturnItemInfoForItem(dict, patchName.c_str(), dict, region.c_str(), ItemType::FACEZONE, nullptr);
}

void ItemInfoFactory::createItemInfoForFileSourceItems
    (
        const string &region,
        const dictionary &fileSourceDict,
        const string &fileSourceItemName
    )
{
    ItemType type = fileSourceDict.lookupOrDefault(idKeys::TYPE_KEY, ItemType(ItemType::FIELD_SAMPLING_FILE_SOURCE));
    createStoreAndReturnItemInfoForItem(
        fileSourceDict,
        fileSourceItemName.c_str(),
        fileSourceDict,
        region.c_str(),
        type.getValue(),
        nullptr
    );
}

void ItemInfoFactory::createItemInfoForInternalBoundary
    (
        const dictionary &dict,
        const string &region,
        const string &patchName
    )
{
    storage_.declareNewItemInfo(
        std::make_unique<InternalBoundaryInfo>(
            patchName,
            region,
            dict,
            meshes_.getMesh(region)));
}

void ItemInfoFactory::createItemInfoForProcessorBoundary
    (
        const string &name,
        const string &region
    )
{
    createStoreAndReturnItemInfoForItem(
        dictionary::null,
        name.c_str(),
        dictionary::null,
        region.c_str(),
        ItemType::PROCESSOR_BOUNDARY,
        nullptr
    );
}

void ItemInfoFactory::createItemInfoForGeometry
    (
        const word& name,
        const dictionary &itemDataDict,
        const SurfaceType &surfaceType
    )
{
    dictionary itemVisualisationDict(name);
    createItemInfoForGeometry(itemVisualisationDict, itemDataDict, surfaceType);
}

void ItemInfoFactory::createItemInfoForBlockMeshGeometry(
    const BlockMeshPatch& patch,
    const BlockMeshBlock& block,
    const List<point>& vertices
)
{
    switch(patch.type.getValue())
    {
        case BlockMeshPatchType::UNKNOWN:
            break;
        case BlockMeshPatchType::WALL:
            storage_.declareNewItemInfo(std::make_unique<BoundarySurfaceInfo>(patch, block, vertices));
            break;
    }
}

void ItemInfoFactory::createItemInfoForGroup
    (
        const dictionary &dict,
        const dictionary &visualisationDict
    )
{
    createStoreAndReturnItemInfoForItem(dict, nullptr, visualisationDict, nullptr, ItemType::GROUP, nullptr);
}




const ItemInfo* ItemInfoFactory::createStoreAndReturnItemInfoForItem(
    const dictionary &itemDict,
    const char* name,
    const dictionary &itemVisualisationDict,
    const char* region,
    int type,
    const char *sceneName
)
{
    const dictionary& dict = itemDict.empty() ? itemVisualisationDict : itemDict;
    string nameString = name ? name : dict.lookupOrDefault<string>(idKeys::NAME_KEY, Utils::getSubdictName(dict));
    word regionString = region ? region : dict.lookupOrDefault<word>(idKeys::REGION_KEY, regionKeys::DEFAULT_REGION_KEY);
    ItemType itemType = type >= 0 ? ItemType(ItemType::Value(type)) : dict.lookup<ItemType>(idKeys::TYPE_KEY);

    Id id(nameString, regionString, itemType);

    if (nameString.empty())
    {
        FatalError << "An item with no name was declared in " << dict.name() << endl << abort(FatalError);
        return nullptr;
    }

    const ItemInfo* cItemInfo;
    if (sceneName && storage_.containsBaseItemInfo(id))
    {
        ItemInfo* itemInfo = storage_.createDeclareAndReturnCopyOfBaseItemInfo(id, sceneName);
        itemInfo->readVisualisationDict(itemVisualisationDict);
        itemInfo->updateObjectData(itemDict);
        cItemInfo = itemInfo;
    }
    else
    {
        std::unique_ptr<ItemInfo> itemInfo = createItemInfoPointerForItem(id, itemVisualisationDict, itemDict);
        if (itemInfo == nullptr)
        {
            itemInfo = std::make_unique<MissingPostDictItemInfo>(id);
        }
        if (sceneName)
        {
            cItemInfo = storage_.declareNewItemInfo(std::move(itemInfo), sceneName);
        }
        else
        {
            cItemInfo = storage_.declareNewItemInfo(std::move(itemInfo));
        }
    }
    return cItemInfo;
}

std::unique_ptr<ItemInfo> ItemInfoFactory::createItemInfoPointerForItem(
    const Id& itemId,
    const dictionary &itemVisualisationDict,
    const dictionary &itemDict
)
{
    switch (itemId.type.getValue())
    {
        case ItemType::PATCH:
            if (regionKeys::VOLUME_MESH_NAME == itemId.name)
            {
                return std::make_unique<VolMeshInfo>(itemId.name, itemId.region, itemVisualisationDict, meshes_.getMesh(itemId.region));
            }
            else
            {
                return std::make_unique<PatchInfo>(itemId.name, itemId.region, itemVisualisationDict, meshes_.getMesh(itemId.region));
            }
        case ItemType::FEATURE_LINE:
        case ItemType::SURFACE:
            return createItemInfoPointerForGeometryItem(
                itemDict,
                itemVisualisationDict,
                itemDict.lookupOrDefault(
                    surfaceKeys::TYPE_KEY,
                    SurfaceType(SurfaceType::UNKNOWN)));
        case ItemType::GROUP:
            return std::make_unique<GroupInfo>(itemDict);
        case ItemType::SLICE:
            return std::make_unique<SliceInfo>(itemDict, itemVisualisationDict, caseFolder_, referenceFrames_);
        case ItemType::CLIP:
            return std::make_unique<ClipInfo>(itemDict, itemVisualisationDict, caseFolder_, referenceFrames_);
        case ItemType::TRANSFORM:
            return std::make_unique<TransformInfo>(itemDict, itemVisualisationDict);
        case ItemType::ISO_SURFACE:
            return std::make_unique<IsoSurfaceInfo>(itemDict, itemVisualisationDict);
        case ItemType::STREAMLINES:
            return std::make_unique<StreamlinesInfo>(itemDict, itemVisualisationDict);
        case ItemType::GLYPH:
            return std::make_unique<GlyphInfo>(itemDict, itemVisualisationDict);
        case ItemType::THRESHOLD:
            return std::make_unique<ThresholdInfo>(itemDict, itemVisualisationDict);
        case ItemType::FIELD_SAMPLING:
            return std::make_unique<FieldSamplingInfo>(itemDict, itemVisualisationDict, caseFolder_);
        case ItemType::FIELD_SAMPLING_FILE_SOURCE:
            return std::make_unique<FieldSamplingFunctionObjectInfo>(itemId.name, itemId.region, caseFolder_, dictionaries_.getFunctionObjectsDict().subDict(itemId.name));
        case ItemType::CELLZONE:
            return std::make_unique<CellZoneInfo>(itemId.name, itemId.region, itemVisualisationDict, meshes_.getMesh(itemId.region));
        case ItemType::FACEZONE:
            return std::make_unique<FaceZoneInfo>(itemId.name, itemId.region, itemVisualisationDict, meshes_.getMesh(itemId.region));
        case ItemType::PROCESSOR_BOUNDARY:
            return std::make_unique<ProcessorBoundaryInfo>(itemId.name, itemId.region, meshes_.getMesh(itemId.region));
        case ItemType::LINE_3D:
            return std::make_unique<Line3DInfo>(itemDict, itemVisualisationDict, caseFolder_);
        case ItemType::MERIDIONAL_GRID:
            return std::make_unique<MeridionalGridInfo>(itemDict, itemVisualisationDict);
        case ItemType::TURBO_SLICE_STREAMWISE:
            return std::make_unique<TurboSliceStreamwiseInfo>(itemDict, itemVisualisationDict);
        case ItemType::TURBO_SLICE_SPANWISE:
            return std::make_unique<TurboSliceSpanwiseInfo>(itemDict, itemVisualisationDict);
        case ItemType::INLET_TO_OUTLET_LINE:
            return std::make_unique<InletToOutletLineInfo>(itemDict, itemVisualisationDict);
        case ItemType::HUB_TO_SHROUD_LINE:
            return std::make_unique<HubToShroudLineInfo>(itemDict, itemVisualisationDict);
        case ItemType::BLADE_LOADING_LINE:
            return std::make_unique<BladeLoadingLineInfo>(itemDict, itemVisualisationDict);
        case ItemType::PVD_FILE_SOURCE:
            return std::make_unique<PvdFileInfo>(itemDict, itemVisualisationDict, itemId.region, caseFolder_);
        case ItemType::IMPORT_DATA:
            return std::make_unique<ImportDataInfo>(itemDict, itemVisualisationDict);
        case ItemType::TOOL_3D_TO_2D:
            return std::make_unique<Tool3Dto2DInfo>(itemDict, itemVisualisationDict);
        case ItemType::SURFACE_SPLITTER_TEST:
            return std::make_unique<SurfaceSplitterInfo>(itemDict, itemVisualisationDict);
        default:
            WarningInFunction << "Item \"" << itemId.name <<
                              "\" has an unknown type - it will not be parsed" << endl;
            return nullptr;
    }
}

std::unique_ptr<ItemInfo> ItemInfoFactory::createItemInfoPointerForGeometryItem(
    const dictionary &itemDict,
    const dictionary &itemVisualisationDict,
    const SurfaceType& surfaceType
)
{
    switch(surfaceType.getValue())
    {
        case SurfaceType::FROM_FILE:
            return std::make_unique<SurfaceFromFileInfo>(itemVisualisationDict, itemDict, caseFolder_);
        case SurfaceType::BOX:
            return std::make_unique<SurfaceBoxInfo>(itemVisualisationDict, itemDict);
        case SurfaceType::SPHERE:
            return std::make_unique<SurfaceSphereInfo>(itemVisualisationDict, itemDict);
        case SurfaceType::CYLINDER:
            return std::make_unique<SurfaceCylinderInfo>(itemVisualisationDict, itemDict);
        case SurfaceType::PLANE:
            return std::make_unique<SurfacePlaneInfo>(itemVisualisationDict, itemDict, meshes_.getDomainBounds());
        case SurfaceType::RING:
            return std::make_unique<SurfaceRingInfo>(itemVisualisationDict, itemDict);
        case SurfaceType::FEATURE_LINE:
            return std::make_unique<FeatureLineInfo>(itemVisualisationDict, caseFolder_);
        case SurfaceType::UNKNOWN:
        default:
            WarningInFunction << "Geometry item " << Id(itemDict.empty() ? itemVisualisationDict : itemDict).name <<
                              " has an unknown surface type - it will not be parsed" << endl;
            return nullptr;
    }
}

void ItemInfoFactory::createItemInfoForGeometry
    (
        const dictionary &itemVisualisationDict,
        const dictionary &itemDataDict,
        const SurfaceType &surfaceType
    )
{
    std::unique_ptr<ItemInfo> itemInfo = createItemInfoPointerForGeometryItem(
        itemDataDict,
        itemVisualisationDict,
        surfaceType
    );
    if (itemInfo == nullptr)
    {
        string dictName = Utils::getSubdictName(itemVisualisationDict);
        storage_.declareNewItemInfo(std::make_unique<MissingGeometryItemInfo>(dictName, ItemType::FEATURE_LINE));
        storage_.declareNewItemInfo(std::make_unique<MissingGeometryItemInfo>(dictName, ItemType::SURFACE));
    }
    else
    {
        storage_.declareNewItemInfo(std::move(itemInfo));
    }
}

// ************************************************************************* //
}