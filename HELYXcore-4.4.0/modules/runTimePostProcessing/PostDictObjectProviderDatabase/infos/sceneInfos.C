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
#include "sceneInfos.H"

#include "postDict/postDictKeys.H"

namespace Foam::functionObjects::runTimeVis
{



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void SceneInfos::readDictionaryAndCreateItemInfos(const dictionary &postDict, ItemInfoFactory &factory,
                                                  const ReferenceFrames& rfs)
{
    scenesInfos_.clear();
    if (postDict.found(sceneKeys::SCENES_DICT_KEY))
    {
        const dictionary &scenesDict = postDict.optionalSubDict(sceneKeys::SCENES_DICT_KEY);

        for (const word &sceneName: scenesDict.toc())
        {
            const dictionary &sceneDict = scenesDict.subDict(sceneName);
            scenesInfos_.emplace_back(sceneDict, factory, rfs);
        }
    }
}

ItemRequirements SceneInfos::getItemRequirements() const
{
    ItemRequirements itemRequirements;
    for (const SceneInfo& sceneInfo : scenesInfos_)
    {
        itemRequirements.mergeFromDownstream(sceneInfo.getItemRequirements());
    }
    return itemRequirements;
}

size_t SceneInfos::getSceneCount() const
{
    return scenesInfos_.size();
}

const SceneInfo &SceneInfos::getSceneInfo(size_t index) const
{
    return scenesInfos_.at(index);
}

const SceneInfo &SceneInfos::getSceneInfo(const std::string &sceneName) const
{
    size_t sceneCount = getSceneCount();
    for (size_t i = 0; i < sceneCount; i++)
    {
        const SceneInfo& sceneInfo = getSceneInfo(i);
        if (sceneInfo.name == sceneName)
        {
            return sceneInfo;
        }
    }
    FatalError << "Tried to access scene " << sceneName << " but it does not exist in the postDict"
               << abort(FatalError);

    return scenesInfos_.at(0);
}


// ************************************************************************* //
}
