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
    Foam::functionObjects::runTimeVis::ColorLegend

Description
    Renders a scalarBar on a runTimeVisualisation scene.

SourceFiles
    colorLegend.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimeVis_colorLegend_H
#define functionObjects_runTimeVis_colorLegend_H

#include "db/dictionary/dictionary.H"
#include "dataStructs/widgets/colourLegendData.H"

#include "vtkSmartPointer.h"
#include "vtkContext2DScalarBarActor.h"
#include "vtkContext2DScalarBarRepresentation.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class vtkDiscretizableColorTransferFunction;
class vtkRenderer;


namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
                          Class scalarBar Declaration
\*---------------------------------------------------------------------------*/

class ColorLegend
{
    const ColourLegendData& legendData_;

    vtkSmartPointer<vtkContext2DScalarBarActor> actor_;

    vtkSmartPointer<vtkContext2DScalarBarRepresentation> representation_;

    void initialiseActor();

public:

    //- Construct from colourLegendData struct
    explicit ColorLegend(const ColourLegendData& legendInfo);

    void updateLookupTable(vtkDiscretizableColorTransferFunction* lut, const scalarMinMax& domainRange);

    void addToRenderer(vtkRenderer* renderer) const;
    void removeFromRenderer(vtkRenderer* renderer) const;

    [[nodiscard]] foamField getName() const {return legendData_.name;};

    [[nodiscard]] bool isVisible() const {return legendData_.visible;};

private:
    [[nodiscard]] bool isVertical() const;
    static string decideAutoRangesFormat(scalar min, scalar max) ;
    static bool shouldUseScientificNotation(scalar range) ;
    void updatePositioning(vtkContext2DScalarBarRepresentation* representation) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
