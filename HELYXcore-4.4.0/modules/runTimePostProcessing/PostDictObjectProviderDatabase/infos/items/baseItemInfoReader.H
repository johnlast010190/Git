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
    (c) 2022-2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::BaseItemInfoReader

Group
    grpGraphicsFunctionObjects

Description
    Reads the base items of the postDict (that don't belong to any scene)

SourceFiles
    baseItemInfoReader.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "db/dictionary/dictionary.H"
#include "primitives/strings/word/word.H"
#include "infos/sceneInfo.H"
#include "infos/items/itemInfoFactory.H"

#include <set>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{
class Dictionaries;
class SurfaceType;
struct FoamMeshInfo;


/*---------------------------------------------------------------------------*\
                      Class BaseItemInfoReader Declaration
\*---------------------------------------------------------------------------*/

class BaseItemInfoReader
{

private:

    BaseItemInfoReader() = default;

public:

    static void readBaseItemsToFactory
        (
            ItemInfoFactory &factory,
            const FoamMeshes &meshes,
            const Dictionaries &dictionaries
        );

private:

    static void createSurfaces(
        ItemInfoFactory &factory,
        const Dictionaries &dictionaries
    );
    static std::set<FoamMeshInfo> findMeshPatches(const FoamMeshes &meshes);
    static void createFileSourcePatches(
        ItemInfoFactory &factory,
        const std::set<FoamMeshInfo>& meshes,
        const dictionary &regionsDict
    );
    static void createProcessorBoundaryPatches(ItemInfoFactory& factory, const FoamMeshes& meshes);
    static void createObjects(ItemInfoFactory& factory, const dictionary& objectsDict, const dictionary& objectsVisualisationDict);
    static void createGroups(ItemInfoFactory& factory, const dictionary& groupsDict, const dictionary& visualisationDict);

    static void createHexMeshDictSurfaces(
        ItemInfoFactory &factory,
        const dictionary &hexMeshDict,
        const dictionary& surfacesVisualisationsDict
    );
    static void createGeometryDictSurfaces(
        ItemInfoFactory &factory,
        const dictionary &geometryDict,
        const dictionary& surfacesVisualisationsDict
    );
    static void createBlockMeshSurfaces(ItemInfoFactory &factory, const dictionary &blockMeshDict);
    static void registerGeometry(
        ItemInfoFactory &factory,
        const dictionary &surfacesVisualisationsDict,
        const dictionary &geometryItemDict,
        const SurfaceType &surfaceType
    );
    static void registerGeometry(
        const word& name,
        ItemInfoFactory &factory,
        const dictionary &surfacesVisualisationsDict,
        const dictionary &geometryItemDict,
        const SurfaceType &surfaceType
    );
    static void createDictSurfaces(
        ItemInfoFactory &factory,
        const Dictionaries &dictionaries
    );
    static void createRegularDictSurfaces(
        ItemInfoFactory &factory,
        const dictionary &surfacesVisualisationsDict,
        const dictionary &surfacesDict
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
