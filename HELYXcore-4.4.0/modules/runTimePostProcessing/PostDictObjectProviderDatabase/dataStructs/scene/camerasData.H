/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : Dev
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
    (c) 2022 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::CamerasData

Description
    Information about the cameras
    Contains the data from the dictionary

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "postDict/postDictKeys.H"
#include "dataStructs/scene/cameraData.H"
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{
class ReferenceFrames;

struct CamerasData
{
    std::vector<CameraData> cameras;

    CamerasData() = default;

    explicit CamerasData(const ReferenceFrames& rfs, const dictionary& sceneDict)
    {
        readDict(rfs, sceneDict);
    }

    void readDict(const ReferenceFrames& rfs, const dictionary& sceneDict)
    {
        cameras.clear();
        if (sceneDict.found(sceneKeys::CAMERAS_DICT_KEY))
        {
            auto camerasDictList = sceneDict.lookup<List<dictionary>>(sceneKeys::CAMERAS_DICT_KEY);
            for (const dictionary &cameraDict: camerasDictList)
            {
                cameras.emplace_back(rfs, cameraDict);
            }
        }
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
