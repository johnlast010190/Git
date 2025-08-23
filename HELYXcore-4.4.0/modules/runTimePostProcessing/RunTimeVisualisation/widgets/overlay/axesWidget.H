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
    Foam::functionObjects::runTimeVis::AxesWidget

Description
    Widget that shows the axes directions (X, Y, Z)

SourceFiles
    axesWidget.C

\*---------------------------------------------------------------------------*/

#pragma once

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "vtkSmartPointer.h"
#include <memory>

class vtkOrientationMarkerWidget;
class vtkRenderWindowInteractor;
class vtkRenderer;

namespace Foam::functionObjects::runTimeVis
{

struct AxisWidgetData;

/*---------------------------------------------------------------------------*\
                        Class AxisWidget Declaration
\*---------------------------------------------------------------------------*/

class AxesWidget
{
private:

    // Private data
    vtkSmartPointer<vtkOrientationMarkerWidget> axesWidget_;

    explicit AxesWidget(const AxisWidgetData& data);

public:
    static std::unique_ptr<AxesWidget> createWidgetIfVisible(const AxisWidgetData& data);


    void addToRenderer(vtkRenderer* overlayRenderer);
    void removeFromRenderer();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
