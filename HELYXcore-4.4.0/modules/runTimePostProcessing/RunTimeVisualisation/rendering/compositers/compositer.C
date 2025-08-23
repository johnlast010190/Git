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

\*---------------------------------------------------------------------------*/

// OpenFOAM includes
#include "compositer.H"

#include "rendering/x3dWriters/serialX3dWriter.H"

#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkSynchronizedRenderWindows.h"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

vtkSmartPointer<vtkSynchronizedRenderWindows> Compositer::createWindowSynchronizer(vtkRenderWindow* window)
{
    vtkSmartPointer<vtkSynchronizedRenderWindows> syncWindows;
    if (controller_)
    {
        syncWindows = vtkSmartPointer<vtkSynchronizedRenderWindows>::New();
        syncWindows->SetRenderWindow(window);
        syncWindows->SetIdentifier(1);
        syncWindows->RenderEventPropagationOff();

        // The GUI does what's in the comment (call render on each process)
        // Probably not a problem
        // false = Call Render() manually on each process - don't use RMI
        syncWindows->SetParallelRendering(true);
        syncWindows->SetParallelController(controller_);
    }
    return syncWindows;
}

void Compositer::renderAndWriteImage
    (
        vtkRenderWindow* window,
        const RendererExtraData& extraData
    )
{
    renderAndWriteImage_(window, extraData);
}

void Compositer::renderAndWriteImage_
    (
        vtkRenderWindow* window,
        const RendererExtraData& extraData
    )
{
    rendererWriter->renderAndWrite(window, extraData);
}

void Compositer::renderAndWriteX3D
(
        vtkRenderWindow* window,
        const fileName& pathToOutputFileWithoutExtension
)
{
    SerialX3dWriter::writeX3D(window, renderer_, pathToOutputFileWithoutExtension);
}

void Compositer::setRendererWriter(RendererWriter* rendererWriter_)
{
    this->rendererWriter.reset(rendererWriter_);
}

} // End namespace



// ************************************************************************* //
