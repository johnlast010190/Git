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
    Foam::functionObjects::runTimeVis::LogoWidget

Description
    Widget that shows the Engys logo

SourceFiles
    logoWidget.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimeVis_logoWidget_H
#define functionObjects_runTimeVis_logoWidget_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "vtkSmartPointer.h"
#include <memory>

class vtkActor2D;
class vtkRenderer;

namespace Foam::functionObjects::runTimeVis
{

struct LogoWidgetData;

/*---------------------------------------------------------------------------*\
                        Class LogoWidget Declaration
\*---------------------------------------------------------------------------*/

class LogoWidget
{
private:

    // Private data

    vtkSmartPointer<vtkActor2D> actor_;
    bool initialised_ = false;

    LogoWidget() = default;
    void initialise(int imageWidth);

public:
    static std::unique_ptr<LogoWidget> createWidgetIfVisible(const LogoWidgetData& data);

    void addToRenderer(vtkRenderer* renderer);
    void removeFromRenderer(vtkRenderer* renderer);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
