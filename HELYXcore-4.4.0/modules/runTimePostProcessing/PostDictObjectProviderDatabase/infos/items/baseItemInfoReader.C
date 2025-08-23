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
    (c) 2023-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include <utility>

#include "baseItemInfoReader.H"

#include "db/dictionary/dictionary.H"
#include "meshes/primitiveShapes/point/point.H"
#include "dictionaries/dictionaries.H"

#include "postDict/postDictKeys.H"
#include "types/surfaceType.H"
#include "types/surfaceFormat.H"

#include "dataStructs/geometry/blockMesh/blockMeshBlock.H"
#include "dataStructs/geometry/blockMesh/blockMeshPatch.H"

#include "engysParallelUtils.h"
#include "engysParallelClass.h"
#include "vtkMultiProcessController.h"

#include "Utils/patchTypes.H"

namespace Foam
{
class BlockMeshPatch
{
public:
    word type_;
    word name_;
    List<List<label>> faces_;
public:
    BlockMeshPatch() = default;

    explicit BlockMeshPatch(Istream &is)
    {
        read(is);
    }

    void read(Istream &is)
    {
        is >> type_;
        is >> name_;
        is >> faces_;

        // Check state of Istream
        is.check(FUNCTION_NAME);
    }
};

Istream &operator>>(Istream &is, BlockMeshPatch &m)
{
    m.read(is);
    return is;
}

class BlockMeshBlock
{
public:
    word type_;
    List<label> points;
    Vector<label> elements;
    word style;
    Vector<label> scale;
public:
    BlockMeshBlock() = default;

    explicit BlockMeshBlock(Istream &is)
    {
        read(is);
    }

    void read(Istream &is)
    {
        is >> type_;
        is >> points;
        is >> elements;
        is >> style;
        is >> scale;

        // Check state of Istream
        is.check(FUNCTION_NAME);
    }
};

Istream &operator>>(Istream &is, BlockMeshBlock &m)
{
    m.read(is);
    return is;
}

namespace functionObjects::runTimeVis
{

struct FoamMeshInfo
{
    enum MeshType {
        PATCH,
        FACE_ZONE,
        CELL_ZONE,
        INTERNAL
    };
    std::string name;
    std::string region;
    MeshType type;

    FoamMeshInfo(word name, string region, MeshType type) : name(std::move(name)), region(std::move(region)), type(type) {};

    bool operator< (const FoamMeshInfo& b) const
    {
        if (this->region != b.region)
        {
            return this->region < b.region;
        }
        else if (this->name != b.name)
        {
            return this->name < b.name;
        }
        else
        {
            return this->type < b.type;
        }
    }

    bool operator==(const FoamMeshInfo& b) const
    {
        return this->region == b.region &&
               this->name == b.name &&
               this->type == b.type;
    }


};

class FoamMeshParallelClass : public engysParallelClass
{
    vtkTypeMacro(FoamMeshParallelClass, engysParallelClass);

    // Description:
    // Construct object
    static FoamMeshParallelClass *New();

public:
    void tradeWithProcess(engysParallelUtils* utils, int remoteProcessId) override
    {
        if (utils->GetController()->GetLocalProcessId() < remoteProcessId)
        {
            receiveFromProcess(utils, remoteProcessId);
            sendToProcess(utils, remoteProcessId);
        }
        else
        {
            sendToProcess(utils, remoteProcessId);
            receiveFromProcess(utils, remoteProcessId);
        }
    }

    void setMeshInfosPointer(std::set<FoamMeshInfo>* infos)
    {
        this->meshInfos = infos;
    }

protected:
    std::set<FoamMeshInfo>* meshInfos;

    void sendToProcess(engysParallelUtils* utils, int remoteProcessId)
    {
        size_t infosToSend = meshInfos->size();
        utils->GetController()->Send(&infosToSend, 1, remoteProcessId, 732);
        for (const FoamMeshInfo& meshInfo : *meshInfos)
        {
            sendMeshInfoToProcess(utils, meshInfo, remoteProcessId);
        }
    }

    void sendMeshInfoToProcess(engysParallelUtils* utils, const FoamMeshInfo& meshInfo, int remoteProcessId)
    {
        utils->sendStringToProcess(meshInfo.region, remoteProcessId);
        utils->sendStringToProcess(meshInfo.name, remoteProcessId);
        int type = meshInfo.type;
        utils->GetController()->Send(&type, 1, remoteProcessId, 614);
    }

    void receiveFromProcess(engysParallelUtils* utils, int remoteProcessId)
    {
        size_t infosToReceive;
        utils->GetController()->Receive(&infosToReceive, 1, remoteProcessId, 732);
        for (size_t i = 0; i < infosToReceive; i++)
        {
            meshInfos->insert(receiveMeshInfoFromProcess(utils, remoteProcessId));
        }
    }

    FoamMeshInfo receiveMeshInfoFromProcess(engysParallelUtils* utils, int remoteProcessId)
    {
        string region = utils->receiveStringFromProcess(remoteProcessId);
        word name = utils->receiveStringFromProcess(remoteProcessId);
        int type;
        utils->GetController()->Receive(&type, 1, remoteProcessId, 614);
        return {name, region, (FoamMeshInfo::MeshType)type};
    }
};
vtkStandardNewMacro(FoamMeshParallelClass)


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void BaseItemInfoReader::readBaseItemsToFactory
    (
        ItemInfoFactory &factory,
        const FoamMeshes &meshes,
        const Dictionaries &dictionaries
    )
{
    // Split dictionary into subdicts
    // (fail if any of these subdicts aren't found)
    const dictionary &visualisationsDict = GET_OPTIONAL_DICTIONARY(dictionaries.getRtppDict(), postDictKeys::VISUALISATIONS_SUBDICT_KEY);
    const dictionary &regionsDict = GET_OPTIONAL_DICTIONARY(visualisationsDict,visualisationKeys::REGIONS_SUBDICT_KEY);

    createSurfaces(factory, dictionaries);

    std::set<FoamMeshInfo> meshInfos = findMeshPatches(meshes);

    vtkMultiProcessController *controller = vtkMultiProcessController::GetGlobalController();
    if (controller && controller->GetNumberOfProcesses() > 1)
    {
        vtkNew<engysParallelUtils> helper;
        helper->SetController(controller);
        vtkNew<FoamMeshParallelClass> parallelClass;
        parallelClass->setMeshInfosPointer(&meshInfos);
        helper->TradeBetweenProcesses(parallelClass);
    }

    createFileSourcePatches(factory, meshInfos, regionsDict);

    createProcessorBoundaryPatches(factory, meshes);

    const dictionary &objectsDict = GET_OPTIONAL_DICTIONARY(dictionaries.getRtppDict(), postDictKeys::OBJECTS_SUBDICT_KEY);
    const dictionary &objectsVisualisationDict = GET_OPTIONAL_DICTIONARY(visualisationsDict, visualisationKeys::OBJECTS_SUBDICT_KEY);
    createObjects(factory, objectsDict, objectsVisualisationDict);

    const dictionary &groupsDict = GET_OPTIONAL_DICTIONARY(dictionaries.getRtppDict(), postDictKeys::GROUPS_SUBDICT_KEY);
    const dictionary &groupsVisualisationDict = GET_OPTIONAL_DICTIONARY(visualisationsDict, visualisationKeys::GROUPS_SUBDICT_KEY);
    createGroups(factory, groupsDict, groupsVisualisationDict);

    factory.createMeshBoundaryGroups();
}

void BaseItemInfoReader::createSurfaces(
    ItemInfoFactory &factory,
    const Dictionaries &dictionaries
)
{
    createDictSurfaces(factory, dictionaries);
    createBlockMeshSurfaces(factory, dictionaries.getBlockMeshDict());
}

void BaseItemInfoReader::createDictSurfaces(
    ItemInfoFactory &factory,
    const Dictionaries &dictionaries
)
{
    if (dictionaries.getRtppDict().isDict(postDictKeys::VISUALISATIONS_SUBDICT_KEY))
    {
        const dictionary& visualisationDict = dictionaries.getRtppDict().subDict(postDictKeys::VISUALISATIONS_SUBDICT_KEY);
        if (visualisationDict.isDict(visualisationKeys::GEOMETRIES_SUBDICT_KEY))
        {
            const dictionary &surfacesVisualisationsDict = visualisationDict.subDict(visualisationKeys::GEOMETRIES_SUBDICT_KEY);
            createGeometryDictSurfaces(factory, dictionaries.getGeometryDict(), surfacesVisualisationsDict);
            createHexMeshDictSurfaces(factory, dictionaries.getHexMeshDict(), surfacesVisualisationsDict);
            return;
        }
    }
    dictionary surfacesVisualisationsDict;
    createGeometryDictSurfaces(factory, dictionaries.getGeometryDict(), surfacesVisualisationsDict);
    createHexMeshDictSurfaces(factory, dictionaries.getHexMeshDict(), surfacesVisualisationsDict);
}

void BaseItemInfoReader::createHexMeshDictSurfaces(
    ItemInfoFactory &factory,
    const dictionary &hexMeshDict,
    const dictionary& surfacesVisualisationsDict
)
{
    if (!hexMeshDict.isDict(surfaceKeys::HEX_MESH_GEOMETRY_SUBDICT_KEY)) return;

    const dictionary &hexMeshGeometryDict = hexMeshDict.subDict(surfaceKeys::HEX_MESH_GEOMETRY_SUBDICT_KEY);
    createRegularDictSurfaces(factory, surfacesVisualisationsDict, hexMeshGeometryDict);

    if (hexMeshDict.found(surfaceKeys::CASTELLATED_MESH_CONTROLS_SUBDICT_KEY))
    {
        const dictionary &hexMeshCastellatedDict = hexMeshDict.subDict(
            surfaceKeys::CASTELLATED_MESH_CONTROLS_SUBDICT_KEY
        );

        if (hexMeshCastellatedDict.found(surfaceKeys::FEATURES_LIST_KEY))
        {
            auto featureDictList = hexMeshCastellatedDict.lookup<List<dictionary>>(
                surfaceKeys::FEATURES_LIST_KEY
            );
            for (const dictionary& featureLineDict : featureDictList)
            {
                const auto nameWithExtension = featureLineDict.lookupOrDefault<string>(surfaceKeys::FILE_KEY, "");
                SurfaceFormat lineFormat = SurfaceFormat::fromNameWithExtension(nameWithExtension);
                word name = nameWithExtension.substr(0, nameWithExtension.length() - lineFormat.getExtension().length());

                registerGeometry(
                    name,
                    factory,
                    surfacesVisualisationsDict,
                    featureLineDict,
                    SurfaceType(SurfaceType::FEATURE_LINE)
                    );
            }
        }
    }
}

void BaseItemInfoReader::createRegularDictSurfaces(
    ItemInfoFactory &factory,
    const dictionary &surfacesVisualisationsDict,
    const dictionary &surfacesDict
)
{
    for (const word& entryName : surfacesDict.toc())
    {
        const dictionary &geometryItemDict = surfacesDict.subDict(entryName);
        auto surfaceType = geometryItemDict.lookupOrDefault(surfaceKeys::TYPE_KEY, SurfaceType(SurfaceType::FROM_FILE));
        if (surfaceType.getValue() == SurfaceType::FROM_FILE)
        {
            word fileName = geometryItemDict.lookupOrDefault(surfaceKeys::FILE_KEY, entryName);
            SurfaceFormat surfaceFormat = SurfaceFormat::fromNameWithExtension(fileName);
            word name = fileName.substr(0, fileName.length() - surfaceFormat.getExtension().length());
            registerGeometry(name, factory, surfacesVisualisationsDict, geometryItemDict, SurfaceType(SurfaceType::FROM_FILE));

            if (geometryItemDict.found(surfaceKeys::SOLIDS_KEY))
            {
                auto solids = geometryItemDict.lookup<List<string>>(surfaceKeys::SOLIDS_KEY);
                for (const string& solid : solids)
                {
                    registerGeometry(name + "_" + solid, factory, surfacesVisualisationsDict, geometryItemDict, SurfaceType(SurfaceType::FROM_FILE));
                }
            }
        }
        else
        {
            registerGeometry(factory, surfacesVisualisationsDict, geometryItemDict, surfaceType);
        }
    }
}

void BaseItemInfoReader::registerGeometry(
    ItemInfoFactory &factory,
    const dictionary &surfacesVisualisationsDict,
    const dictionary &geometryItemDict,
    const SurfaceType &surfaceType
)
{
    const word geometryName = geometryItemDict.lookupOrDefault<string>(
        surfaceKeys::NAME_KEY,
        Utils::getSubdictName(geometryItemDict));
    registerGeometry(geometryName, factory, surfacesVisualisationsDict, geometryItemDict, surfaceType);
}

void BaseItemInfoReader::registerGeometry(
    const word& name,
    ItemInfoFactory &factory,
    const dictionary &surfacesVisualisationsDict,
    const dictionary &geometryItemDict,
    const SurfaceType &surfaceType
)
{
    bool foundVisualisationEntry = false;
    if (!surfacesVisualisationsDict.empty())
    {
        for (const word &surfaceName: surfacesVisualisationsDict.toc())
        {
            const dictionary &itemVisualisationDict = surfacesVisualisationsDict.subDict(surfaceName);

            if (name == surfaceName)
            {
                factory.createItemInfoForGeometry(itemVisualisationDict, geometryItemDict, surfaceType);
                foundVisualisationEntry = true;
            }

            word hexSurfaceName = itemVisualisationDict.lookupOrDefault(surfaceKeys::PARENT_NAME_KEY, surfaceName);
            if (name == hexSurfaceName)
            {
                factory.createItemInfoForGeometry(itemVisualisationDict, geometryItemDict, surfaceType);
            }
        }
    }

    if (!foundVisualisationEntry)
    {
        factory.createItemInfoForGeometry(name, geometryItemDict, surfaceType);
    }
}

void BaseItemInfoReader::createGeometryDictSurfaces(
    ItemInfoFactory &factory,
    const dictionary &geometryDict,
    const dictionary& surfacesVisualisationsDict
)
{
    if (geometryDict.empty()) return;

    const dictionary &surfacesDict = geometryDict.subDict(geometryDictKeys::SURFACES_SUBDICT_KEY);
    createRegularDictSurfaces(factory, surfacesVisualisationsDict, surfacesDict);

    if (geometryDict.found(geometryDictKeys::LINES_SUBDICT_KEY))
    {
        const dictionary &linesDict = geometryDict.subDict(geometryDictKeys::LINES_SUBDICT_KEY);
        for (const word& nameWithExtension : linesDict.toc())
        {
            SurfaceFormat lineFormat = SurfaceFormat::fromNameWithExtension(nameWithExtension);
            word name = nameWithExtension.substr(0, nameWithExtension.length() - lineFormat.getExtension().length());
            const dictionary &geometryItemDict = linesDict.subDict(nameWithExtension);
            registerGeometry(
                name,
                factory,
                surfacesVisualisationsDict,
                geometryItemDict,
                SurfaceType(SurfaceType::FEATURE_LINE)
            );
        }
    }
}

void BaseItemInfoReader::createBlockMeshSurfaces(ItemInfoFactory &factory, const dictionary &blockMeshDict)
{
    if (blockMeshDict.found(blockMeshDictKeys::VERTICES_KEY) &&
        blockMeshDict.found(blockMeshDictKeys::BLOCKS_KEY) &&
        blockMeshDict.found(blockMeshDictKeys::PATCHES_KEY))
    {
        auto vertices = blockMeshDict.lookup<List<vector>>(blockMeshDictKeys::VERTICES_KEY);
        auto blocks = blockMeshDict.lookup<List<BlockMeshBlock>>(blockMeshDictKeys::BLOCKS_KEY);
        auto patches = blockMeshDict.lookup<List<BlockMeshPatch>>(blockMeshDictKeys::PATCHES_KEY);

        for (const BlockMeshPatch &patch: patches)
        {
            factory.createItemInfoForBlockMeshGeometry(patch, blocks[0], vertices);
        }
    }
}

std::set<FoamMeshInfo> BaseItemInfoReader::findMeshPatches(const FoamMeshes &meshes)
{
    std::set<FoamMeshInfo> foamMeshInfos;
    for (const string &region: meshes.getMeshes().toc())
    {
        const fvMesh &mesh = meshes.getMesh(region).getMesh();

        foamMeshInfos.emplace(regionKeys::VOLUME_MESH_NAME, region, FoamMeshInfo::PATCH);

        for (const fvPatch& patch : mesh.boundary())
        {
            if (PatchTypes::isInternalPatch(patch))
            {
                foamMeshInfos.emplace(patch.name(), region, FoamMeshInfo::INTERNAL);
            }
            else if (PatchTypes::isProcessPatch(patch) ||
                     PatchTypes::isNCCPatch(patch))
            {
                // do nothing
            }
            else
            {
                foamMeshInfos.emplace(patch.name(), region, FoamMeshInfo::PATCH);
            }
        }

        for (const faceZone& face : mesh.faceZones())
        {
            foamMeshInfos.emplace(face.name(), region, FoamMeshInfo::FACE_ZONE);
        }

        for (const cellZone& cell : mesh.cellZones())
        {
            foamMeshInfos.emplace(cell.name(), region, FoamMeshInfo::CELL_ZONE);
        }
    }
    return foamMeshInfos;
}

void BaseItemInfoReader::createFileSourcePatches(
    ItemInfoFactory &factory,
    const std::set<FoamMeshInfo>& meshes,
    const dictionary &regionsDict
)
{
    for (const word &region : regionsDict.toc())
    {
        const dictionary &regionDict = regionsDict.subDict(region);
        for (const word &patchName: regionDict.toc())
        {
            const dictionary &patchDict = regionDict.subDict(patchName);
            if (regionKeys::FILE_SOURCES_KEY == patchName)
            {
                for (const word &fileSourceName: patchDict.toc())
                {
                    const dictionary &fileSourceDict = patchDict.subDict(fileSourceName);
                    factory.createItemInfoForFileSourceItems(
                        region,
                        fileSourceDict,
                        fileSourceName
                    );
                }
            }
        }
    }

    for (const FoamMeshInfo& mesh : meshes)
    {
        const dictionary &regionDict = GET_OPTIONAL_DICTIONARY(regionsDict, mesh.region);
        const dictionary &patchDict = GET_OPTIONAL_DICTIONARY(regionDict, mesh.name);

        switch(mesh.type)
        {
            case FoamMeshInfo::PATCH:
                factory.createItemInfoForPatch(patchDict, mesh.region,mesh.name);
                break;
            case FoamMeshInfo::FACE_ZONE:
                factory.createItemInfoForFaceZone(patchDict, mesh.region,mesh.name);
                break;
            case FoamMeshInfo::CELL_ZONE:
                factory.createItemInfoForCellZone(patchDict, mesh.region,mesh.name);
                break;
            case FoamMeshInfo::INTERNAL:
                factory.createItemInfoForInternalBoundary(patchDict, mesh.region,mesh.name);
                break;
        }
    }
}

void BaseItemInfoReader::createProcessorBoundaryPatches(ItemInfoFactory &factory, const FoamMeshes &meshes)
{
    for (const string &meshName: meshes.getMeshes().toc())
    {
        const fvMesh &mesh = meshes.getMesh(meshName).getMesh();
        for (const auto & boundary : mesh.boundary())
        {
            if (PatchTypes::isProcessPatch(boundary))
            {
                factory.createItemInfoForProcessorBoundary(boundary.name(), meshName);
            }
        }
    }
}

void BaseItemInfoReader::createObjects
    (
        ItemInfoFactory &factory,
        const dictionary &objectsDict,
        const dictionary &objectsVisualisationDict
    )
{
    for (const word &objectName: objectsDict.toc())
    {
        const dictionary &objectDict = objectsDict.subDict(objectName);
        const dictionary &objectVisDict = objectsVisualisationDict.isDict(objectName) ?
            objectsVisualisationDict.subDict(objectName) : dictionary::null;

        factory.createItemInfoForObject(objectDict, objectVisDict);
    }
}


void BaseItemInfoReader::createGroups(ItemInfoFactory &factory, const dictionary &groupsDict, const dictionary& visualisationDict)
{
    for (const word &groupName: groupsDict.toc())
    {
        const dictionary &groupDict = groupsDict.subDict(groupName);
        factory.createItemInfoForGroup(groupDict, GET_OPTIONAL_DICTIONARY(visualisationDict, groupName));
    }
}

} // End namespace functionObjects
} // End namespace Foam

// ************************************************************************* //
