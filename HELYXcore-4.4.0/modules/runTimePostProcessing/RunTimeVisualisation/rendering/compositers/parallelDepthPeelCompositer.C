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

\*---------------------------------------------------------------------------*/

#include "compileOptions.H"
#if TRANSPARENT_COMPOSITER_TO_USE == PARALLEL_DEPTH_PEEL_COMPOSITER

// OpenFOAM includes
#include "parallelDepthPeelCompositer.H"

#include "rendering/rtppActor.H"
#include "rendering/pngRenderers/parallelPngWriter.H"

// vtk includes
#include "vtkCompositedSynchronizedRenderers.h"
#include "vtkRenderWindow.h"
#include "vtkParallelOpenGLRenderer.h"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

    ParallelDepthPeelCompositer::ParallelDepthPeelCompositer(vtkMultiProcessController* controller, const RenderInfo& renderInfo)
:
    Compositer(controller, new ParallelPngWriter(renderInfo))
{
    vtkSmartPointer<vtkParallelOpenGLRenderer> compositerRenderer = vtkSmartPointer<vtkParallelOpenGLRenderer>::New();
    compositerRenderer->SetController(controller);
    renderer_ = compositerRenderer;
    synchronizer_ = vtkSmartPointer<vtkSynchronizedRenderers>::New();
    synchronizer_->ParallelRenderingOn();
    synchronizer_->SetParallelController(controller_);
}

void ParallelDepthPeelCompositer::renderAndWriteImage
(
        vtkRenderWindow* window,
        const RendererExtraData& extraData
)
{
    synchronizer_->SetRenderer(renderer_);

    vtkSmartPointer<vtkSynchronizedRenderWindows> windowSync = createWindowSynchronizer(window);
    renderAndWriteImage_(window, extraData);
}

} // End namespace

// ************************************************************************* //
#endif