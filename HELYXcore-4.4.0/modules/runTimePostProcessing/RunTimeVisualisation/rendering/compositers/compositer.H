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
    Foam::functionObjects::runTimeVis::Compositer

Description
    base class for the compositers

SourceFiles
    <compositer.C>

\*---------------------------------------------------------------------------*/

#pragma once

#include "primitives/strings/fileName/fileName.H"
#include "rendering/pngRenderers/pngWriter.H"

#include "vtkSmartPointer.h"
#include <memory>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward declarations
class vtkMultiProcessController;
class vtkRenderer;
class vtkRenderWindow;
class vtkSynchronizedRenderWindows;
class vtkActor;

namespace Foam::functionObjects::runTimeVis
{

class rtppActor;
class RendererWriter;

/*---------------------------------------------------------------------------*\
                         Class Compositer Declaration
\*---------------------------------------------------------------------------*/

class Compositer
{
protected:
    vtkMultiProcessController* controller_;
    vtkSmartPointer<vtkRenderer> renderer_;

public:
    virtual void registerActorPolyData(rtppActor& actor, const char* name) {};
    virtual void registerWidgetActorPolyData(vtkActor* actor) {};

    virtual void redistributeActorPolyData(rtppActor& actor) {};
    virtual void redistributeWidgetActorPolyData(vtkActor* actor) {};

    vtkSmartPointer<vtkRenderer> get3DRenderer() {return renderer_;};

    virtual void renderAndWriteImage
    (
            vtkRenderWindow* window,
            const RendererExtraData& extraData
    );

    virtual void renderAndWriteX3D
    (
            vtkRenderWindow* window,
            const fileName& pathToOutputFileWithoutExtension
    );

    virtual void resetForNewIteration() {};

    virtual ~Compositer() = default;

private:
    std::unique_ptr<RendererWriter> rendererWriter;

protected:
    Compositer(vtkMultiProcessController* controller, RendererWriter* rendererWriter)
        : controller_(controller), rendererWriter(rendererWriter)
    {};
    explicit Compositer(vtkMultiProcessController* controller)
        : controller_(controller)
    {};
    void setRendererWriter(RendererWriter* rendererWriter);

    vtkSmartPointer<vtkSynchronizedRenderWindows> createWindowSynchronizer(vtkRenderWindow* window);

    void renderAndWriteImage_(
        vtkRenderWindow *window,
        const RendererExtraData& extraData
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
