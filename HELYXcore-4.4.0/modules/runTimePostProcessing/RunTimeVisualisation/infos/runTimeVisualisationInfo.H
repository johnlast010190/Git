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
    Foam::functionObjects::runTimeVis::RunTimeVisualisationInfo

Description
    Contains data that is shared by all scenes, and is thus stored in the runTimeVisualisation

SourceFiles
    <runTimeVisualisationInfo.C>

\*---------------------------------------------------------------------------*/

#pragma once


#include "infos/renderInfo.H"
#include "infos/pvdInfo.H"
#include "colourMaps/colourMaps.H"
#include "storage/foamMeshes.H"
#include "storage/itemStorage.H"
#include "colourMaps/colourMaps.H"
#include "infos/items/itemInfoFactory.H"
#include "infos/colourLookupTablesInfo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class RunTimeVisualisationInfo
{
private:

    //- Function name (i.e. the key of the top-level RTPP dictionary)
    // Somehow this goes out of scope sometimes, so take a deep copy
    const word functionName_;

    fileName outputFolder_;
    ColourMaps colourMaps_;
    RenderInfo renderInfo_;
    PvdInfo pvdInfo_;
    ColourLookupTablesInfo defaultColourLookupTablesInfo_;
    List<word> activeScenes_;

public:
    RunTimeVisualisationInfo(
        const word& name,
        const Time& runTime,
        const dictionary& foDict,
        const dictionary& rtppDict
    );

    [[nodiscard]] const word& getFunctionName() const {return functionName_;};
    [[nodiscard]] const RenderInfo& getRenderInfo() const {return renderInfo_;};
    [[nodiscard]] const PvdInfo& getPvdInfo() const {return pvdInfo_;};
    [[nodiscard]] const ColourMaps& getColourMaps() const {return colourMaps_;};
    [[nodiscard]] const fileName& getOutputFolder() const {return outputFolder_;};
    [[nodiscard]] const List<word>& getActiveScenes() const {return activeScenes_;};

    [[nodiscard]] const ColourLookupTablesInfo& getDefaultColourLookupTablesInfo() const {return defaultColourLookupTablesInfo_;};
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
