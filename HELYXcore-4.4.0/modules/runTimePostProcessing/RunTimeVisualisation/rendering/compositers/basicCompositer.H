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
    (c) 2020-2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::BasicCompositer

Description
    Class to synchronize and composite parallel renders.
    Doesn't work with transparency. Does not redistribute data.

SourceFiles
    basicCompositer.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "rendering/compositers/compositer.H"

#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward declarations
class vtkRenderer;
class vtkSynchronizedRenderers;

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
                       Class BasicCompositer Declaration
\*---------------------------------------------------------------------------*/

class BasicCompositer : public Compositer
{
private:
    vtkSmartPointer<vtkSynchronizedRenderers> synchronizer_;

public:

    explicit BasicCompositer(vtkMultiProcessController* controller, const RenderInfo& renderInfo);

    void renderAndWriteImage
    (
            vtkRenderWindow* window,
            const RendererExtraData& extraData
    ) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
