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
    (c) 2024-2025 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::VectorWidget

Description
    Widget that shows an arrow and/or text to represent a constant vector field,
    namely the frame acceleration, though it could be changed to show any vector

SourceFiles
    vectorWidget.C

\*---------------------------------------------------------------------------*/

#pragma once

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "db/Time/Time.H"

#include "vtkSmartPointer.h"
#include "types/vectorDataType.H"
#include <memory>

class engysVectorWidget;
class vtkRenderer;

namespace Foam::functionObjects::runTimeVis
{

struct VectorWidgetData;

/*---------------------------------------------------------------------------*\
                        Class VectorWidget Declaration
\*---------------------------------------------------------------------------*/

class VectorWidget
{
private:

    // Private data
    vtkSmartPointer<engysVectorWidget> vectorWidget_;

    VectorDataType dataSourceType_;
    bool referenceVectorSet_ = false;

    explicit VectorWidget(const VectorWidgetData& data);

public:
    static std::unique_ptr<VectorWidget> createWidgetIfVisible(const VectorWidgetData& data);

    void update(const Time& runTime);

    void addToRenderer(vtkRenderer* overlayRenderer);
    void removeFromRenderer();

    void setReferenceVector(const dictionary &vectorDict);

    void setVectorToShow(const dictionary &vectorDict);

    const dictionary& getSolverDict(const dictionary& resultsDict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
