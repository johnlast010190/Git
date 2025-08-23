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
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::itemInfoFactory

Description
    Class made to return an appropriate itemInfo object for the type

SourceFiles
    itemInfoFactory.C

\*---------------------------------------------------------------------------*/

#pragma once

#include <unordered_map>
#include <utility>
#include "infos/items/itemInfo.H"
#include "storage/foamMeshes.H"
#include "storage/itemStorage.H"
#include "storage/referenceFrames.H"
#include "types/surfaceType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct BlockMeshPatch;
struct BlockMeshBlock;

class ItemInfoFactory
{
private:

    string caseFolder_;
    FoamMeshes &meshes_;
    ItemStorage &storage_;
    const ReferenceFrames& referenceFrames_;
    const Dictionaries& dictionaries_;

public:

    // Constructors

    ItemInfoFactory(
        string caseFolder,
        FoamMeshes &meshes,
        ItemStorage &storage,
        const ReferenceFrames &referenceFrames,
        const Dictionaries &dictionaries
    ) :
        caseFolder_(std::move(caseFolder)),
        meshes_(meshes),
        storage_(storage),
        referenceFrames_(referenceFrames),
        dictionaries_(dictionaries)
    {};

    // Member Functions
    const ItemInfo* createAndReturnItemInfoForSceneItem(const dictionary &dict, const word& sceneName);
    static const dictionary& getSceneObjectItemDict(const dictionary &dict);
    void createItemInfoForObject(const dictionary &dict, const dictionary &objectVisDict);
    void createItemInfoForPatch(const dictionary &dict, const string &region, const string &patchName);
    void createItemInfoForCellZone(const dictionary &dict, const string &region, const string &patchName);
    void createItemInfoForFaceZone(const dictionary &dict, const string &region, const string &patchName);
    void createItemInfoForFileSourceItems(
        const string &region,
        const dictionary &fileSourceVisDict,
        const string &fileSourceItemName
    );

    void createItemInfoForInternalBoundary(const dictionary &dict, const string &region, const string &patchName);
    void createItemInfoForProcessorBoundary(const string &name, const string &region);
    void createItemInfoForGeometry(
        const word& name,
        const dictionary &itemDataDict,
        const SurfaceType &surfaceType
    );
    void createItemInfoForGeometry(
        const dictionary &itemVisualisationDict,
        const dictionary &itemDataDict,
        const SurfaceType &surfaceType
    );
    void createItemInfoForBlockMeshGeometry(
        const BlockMeshPatch& patch,
        const BlockMeshBlock& block,
        const List<point>& vertices
    );
    void createItemInfoForGroup(const dictionary &dict, const dictionary &visualisationDict);

    void createMeshBoundaryGroups()
    { storage_.createMeshBoundaryGroups(); };

private:
    const ItemInfo* createStoreAndReturnItemInfoForItem(
        const dictionary &itemDict,
        const char* name,
        const dictionary &itemVisualisationDict,
        const char* region,
        int type,
        const char *sceneName
    );

    std::unique_ptr<ItemInfo> createItemInfoPointerForItem(
        const Id &itemId,
        const dictionary &itemVisualisationDict,
        const dictionary &itemDict
    );
    std::unique_ptr<ItemInfo> createItemInfoPointerForGeometryItem(
        const dictionary &itemDict,
        const dictionary &itemVisualisationDict,
        const SurfaceType& surfaceType
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
