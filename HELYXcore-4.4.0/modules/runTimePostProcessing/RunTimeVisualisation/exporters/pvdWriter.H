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
    Foam::functionObjects::runTimeVis::pvdWriter

Description
    A class responsible for exporting the contents of a scene to a collection of pvd files

SourceFiles
    pvdWriter.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "infos/sceneInfo.H"
#include "colourMaps/colourMaps.H"
#include "colourLookupTable/colourLookupTableProvider.H"
#include "rendering/renderManager.H"
#include "scene/items.H"
#include <memory>

class vtkMultiProcessController;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
class engysPVDWriter;
namespace Foam
{
class Time;
}

namespace Foam::functionObjects::runTimeVis
{

class RunTimeVisualisationInfo;
class PostDictObjectProviderDatabase;

/*---------------------------------------------------------------------------*\
                          Class PVDWriter Declaration
\*---------------------------------------------------------------------------*/

class PVDWriter
{
private:

    vtkSmartPointer<engysPVDWriter> pvdWriter_;
    Items& items_;
    FoamFields allowableFields_;
    bool writePointData_;

public:

    // Constructors

    PVDWriter
        (
            const std::string&  name,
            const RunTimeVisualisationInfo& rtppInfo,
            const SceneInfo& sceneInfo,
            vtkMultiProcessController* controller,
            Items &items
        );

    // Member Functions

    void exportIfRequired(const Time &currentTime);
    void removeUnwantedFields(vtkDataSetAttributes *data, std::set<std::string> wantedFields);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

