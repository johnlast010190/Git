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
    (c) 2015-2019 OpenCFD Ltd.
    (c) 2022-2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::Scene

Description
    A class that contains everything necessary to render and export a scene

SourceFiles
    scene.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimeVis_scene_H
#define functionObjects_runTimeVis_scene_H

#include "infos/sceneInfo.H"
#include "colourMaps/colourMaps.H"
#include "colourLookupTable/colourLookupTableProvider.H"
#include "rendering/renderManager.H"
#include "scene/items.H"
#include "exporters/pvdWriter.H"
#include <memory>

class vtkMultiProcessController;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
class Time;
}

namespace Foam::functionObjects::runTimeVis
{

class RunTimeVisualisationInfo;
class PostDictObjectProviderDatabase;

/*---------------------------------------------------------------------------*\
                          Class Scene Declaration
\*---------------------------------------------------------------------------*/

class Scene
{
private:

    const std::string name_;
    const SceneInfo &sceneInfo_;

    const PostDictObjectProviderDatabase& database_;
    const RunTimeVisualisationInfo& rtppInfo_;

    std::unique_ptr<Compositer> compositer_;

    ColourLookupTableProvider lookupTableProvider_;

    Items items_;

    const fileName outputDirectory_;

    renderManager renderManager_;
    PVDWriter pvdWriter_;

public:

    // Constructors

    Scene
        (
            std::string sceneName,
            const PostDictObjectProviderDatabase &database,
            const RunTimeVisualisationInfo& rtppInfo
        );

    // Member Functions

    void renderAndExport
        (
            const Time &currentTime,
            const Time &runTime
        );

    std::string getName()
    { return name_; };

private:
    static std::unique_ptr<Compositer> createCompositer(
        vtkMultiProcessController *controller,
        const RenderInfo &renderInfo,
        const HashTable<const ItemInfo *, Id, IdHasher> &itemInfos
    );

    string getOutputFileNameWithoutExtension(const Time &currentTime, const string& cameraName);
    string getOutputFileNameWithoutExtension(const Time &currentTime);
    static string getCameraName(const CameraData& camera, int cameraIndex, size_t totalCameras);
    static string getTimeName(const Time &currentTime);
    void updateCamera(const CameraData &camera, const Time &currentTime, const Time &runTime);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

#endif
