/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : Dev
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
    Foam::functionObjects::runTimeVis::EDFCompositer

Description
    Class to synchronize and composite parallel renders and write an Engys Dynamic Fields file.
    Works with transparency, by applying the depth peeling pass in parallel with all processes.
    Does not redistribute vtk data.

SourceFiles
    EDFCompositer.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "compileOptions.H"

#include "rendering/compositers/compositer.H"

#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward declarations
class vtkRenderer;
class vtkSynchronizedRenderers;
class engysEDFPropsToRender;

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
                   Class EDFCompositer Declaration
\*---------------------------------------------------------------------------*/

class EDFCompositer : public Compositer
{
private:
    vtkSmartPointer<vtkSynchronizedRenderers> synchronizer_;
    vtkSmartPointer<engysEDFPropsToRender> propsToRender_;
    bool useFieldsShownByProps;

public:

    explicit EDFCompositer(vtkMultiProcessController* controller, const RenderInfo& renderInfo);

    void registerActorPolyData(rtppActor& actor, const char* name) override;
    void renderAndWriteImage
    (
            vtkRenderWindow* window,
            const RendererExtraData& extraData
    ) override;

    void resetForNewIteration() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
