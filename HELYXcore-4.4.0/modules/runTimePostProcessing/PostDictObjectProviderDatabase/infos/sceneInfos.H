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
    Foam::functionObjects::runTimeVis::SceneInfos

Description
    Information about the scenes
    Contains the data from the dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "baseClasses/itemRequirements.H"
#include "infos/sceneInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class SceneInfos
{
private:
    std::vector<SceneInfo> scenesInfos_;

public:

    SceneInfos() = default;

    void readDictionaryAndCreateItemInfos(const dictionary& postDict, ItemInfoFactory& factory,
                                          const ReferenceFrames& rfs);

    [[nodiscard]] ItemRequirements getItemRequirements() const;

    [[nodiscard]] size_t getSceneCount() const;
    [[nodiscard]] const SceneInfo &getSceneInfo(size_t index) const ;
    [[nodiscard]] const SceneInfo &getSceneInfo(const std::string& sceneName) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
