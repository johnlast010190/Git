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
    (c) 2020-2021 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::DummyCompositer

Description
    To be used with serial cases, that don't need a compositer or synchronizer

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "rendering/compositers/compositer.H"
#include "rendering/pngRenderers/serialPngWriter.H"

#include "vtkRenderWindow.h"
#include "vtkRenderer.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward declarations

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
                      Class DummyCompositer Declaration
\*---------------------------------------------------------------------------*/

class DummyCompositer : public Compositer
{
public:

    explicit DummyCompositer(const RenderInfo& renderInfo)
        : Compositer(nullptr, new SerialPngWriter(renderInfo))
    {
        renderer_ = vtkSmartPointer<vtkRenderer>::New();
    };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
