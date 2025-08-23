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
    (c) 2019 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::ColorLegends

Description
    Renders the scalarBars on a runTimeVisualisation scene.

SourceFiles
    colorLegends.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "dataStructs/widgets/colourLegendsData.H"
#include "widgets/overlay/colorLegend.H"
#include "storage/externalFields.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class vtkRenderer;

namespace Foam::functionObjects::runTimeVis
{

class ColourLookupTableProvider;
struct ColourMaps;
class FoamMeshes;

/*---------------------------------------------------------------------------*\
                          Class ColorLegends Declaration
\*---------------------------------------------------------------------------*/

class ColorLegends
{
private:
    HashTable<ColorLegend, word> colourLegendScalarBars_;

public:

    static std::unique_ptr<ColorLegends> createWidgetIfVisible
    (
            const ColourLegendsData& data
    );

    void update(
            ColourLookupTableProvider& colourLutProvider,
            const ColourMaps& colourMaps,
            const FoamMeshes& meshes,
            const ExternalFields& externalFields
    );

    void addToRenderer(vtkRenderer* renderer);
    void removeFromRenderer(vtkRenderer* renderer);

private:
    //- Construct from colourLegendData struct
    explicit ColorLegends(const ColourLegendsData& data);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
