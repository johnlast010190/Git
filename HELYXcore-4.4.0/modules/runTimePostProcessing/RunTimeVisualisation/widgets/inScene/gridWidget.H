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
    Foam::functionObjects::runTimeVis::GridWidget

Description
    Widget that shows the grid with the coordinates along the bounds of the domain

SourceFiles
    gridWidget.C

\*---------------------------------------------------------------------------*/

#pragma once

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dataStructs/widgets/gridWidgetData.H"
#include "rendering/compositers/compositer.H"

#include "vtkSmartPointer.h"

class vtkGridAxes3DActor;
class vtkRenderer;

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
                        Class GridWidget Declaration
\*---------------------------------------------------------------------------*/

class GridWidget
{
private:

    const static label FONT_SIZE = 12;

    enum Axes{
        X,
        Y,
        Z
    };

    // Private data

    vtkSmartPointer<vtkGridAxes3DActor> axesActor_;
    GridWidgetData data_;

public:

    // Constructors

        explicit GridWidget(GridWidgetData  data);

        [[nodiscard]] bool isVisible() const {return data_.showAxesGrid;};

        void update(vtkRenderer* renderer);

        void registerPolyDataToCompositer(vtkRenderer* renderer, Compositer* compositer);
        void redistributePolyDataWithCompositer(Compositer* compositer);

        void addToRenderer(vtkRenderer* renderer);
        void removeFromRenderer(vtkRenderer* renderer);

        [[nodiscard]] bool mustUpdateBoundsWithSceneBounds() const;

private:
        void removeLabels();

        void configureCustomLabels(scalar bounds[6]);
        void configureLabelForAxis(const scalar bounds[6], bool useLabelDistance, scalar labelDistance, label axis);

        void configureGridColor();
        void configureLabelColor();
        void configureAxesTitles();
        void configureFaces();
        [[nodiscard]] unsigned int getFaceMaskFromData() const;
        void configureFonts();
        void configureCustomBounds();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
