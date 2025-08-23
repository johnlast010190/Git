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
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::WidgetsInScene

Description
    A class that holds all widgets that have a 3D actor that is added to the scene

SourceFiles
    widgetsInScene.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "widgets/inScene/gridWidget.H"
#include "widgets/inScene/referenceFrameWidget.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct WidgetsInfo;
class ReferenceFrames;

/*---------------------------------------------------------------------------*\
                        Class WidgetsInScene Declaration
\*---------------------------------------------------------------------------*/

class WidgetsInScene
{
private:
    vtkRenderer* baseRenderer_;
    Compositer* compositer_;

    GridWidget gridWidget_;
    std::vector<ReferenceFrameWidget> referenceFrameWidgets_;

public:
    WidgetsInScene(Compositer* compositer, const WidgetsInfo& widgetsInfo, const ReferenceFrames& referenceFrames);

    void setRenderer(vtkRenderer* baseRenderer) {baseRenderer_ = baseRenderer;};

    void addVisibleWidgetsToRenderer();

    void updateVisibleWidgets();

    void registerDataToCompositer();
    void redistributeDataWithCompositer();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
