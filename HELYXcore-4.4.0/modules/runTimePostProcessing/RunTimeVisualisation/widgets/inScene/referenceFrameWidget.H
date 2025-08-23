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
    Foam::functionObjects::runTimeVis::ReferenceFrameWidget

Description
    Widget that shows a reference frame

SourceFiles
    referenceFrameWidget.C

\*---------------------------------------------------------------------------*/

#pragma once

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dataStructs/referenceFrame/referenceFrame.H"
#include "rendering/compositers/compositer.H"

#include "vtkSmartPointer.h"

class engysCoordAxesActor;
class vtkRenderer;

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
                     Class ReferenceFrameWidget Declaration
\*---------------------------------------------------------------------------*/

class ReferenceFrameWidget
{
public:
    // Constructors

    explicit ReferenceFrameWidget(const ReferenceFrame& referenceFrame);

    void update();

    void registerPolyDataToCompositer(Compositer *compositer) const;
    void redistributePolyDataWithCompositer(Compositer *compositer) const;

    void addToRenderer(vtkRenderer *renderer) const;
    void removeFromRenderer(vtkRenderer *renderer) const;

private:

    vtkSmartPointer<engysCoordAxesActor> actor_;
    const ReferenceFrame& referenceFrame_;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
