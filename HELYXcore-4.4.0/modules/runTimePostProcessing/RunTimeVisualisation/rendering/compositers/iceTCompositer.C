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

#include "compileOptions.H"
#if TRANSPARENT_COMPOSITER_TO_USE == ICE_T_COMPOSITER

// OpenFOAM includes
#include "iceTCompositer.H"

#include "rendering/rtppActor.H"
#include "rendering/pngRenderers/parallelPngWriter.H"

// vtk includes
#include "vtkIceTSynchronizedRenderers.h"
#include "vtkPKdTree.h"
#include "vtkKdTreeManager.h"
#include "vtkPartitionOrderingInterface.h"
#include "vtkMapper.h"
#include "vtkOrderedCompositeDistributor.h"
#include "vtkAlgorithmOutput.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderWindow.h"
#include "vtkActor.h"
#include "vtkPolyDataNormals.h"
#include "vtkRenderer.h"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

IceTCompositer::IceTCompositer(vtkMultiProcessController* controller, const RenderInfo& renderInfo)
:
    Compositer(controller, new ParallelPngWriter(renderInfo))
{
    renderer_ = vtkSmartPointer<vtkRenderer>::New();

    synchronizer_ = vtkSmartPointer<vtkIceTSynchronizedRenderers>::New();
    synchronizer_->ParallelRenderingOn();
    synchronizer_->SetUseOrderedCompositing(true);
    synchronizer_->SetParallelController(controller_);

    treeManager_ = vtkSmartPointer<vtkKdTreeManager>::New();
    pKdTree_ = nullptr;

    redistributor_ = vtkSmartPointer<vtkOrderedCompositeDistributor>::New();
    redistributor_->SetController(controller_);
    redistributor_->SetPassThrough(false);
    redistributor_->SetBoundaryMode(vtkOrderedCompositeDistributor::SPLIT_BOUNDARY_CELLS);
}

void IceTCompositer::resetForNewIteration()
{
    treeManager_ = vtkSmartPointer<vtkKdTreeManager>::New();
    pKdTree_ = nullptr;
    redistributor_->SetPKdTree(nullptr);
}

void IceTCompositer::registerActorPolyData(rtppActor& actor, const char* name)
{
    treeManager_->AddDataObject(actor.getProcessedInputData());
}

void IceTCompositer::registerWidgetActorPolyData(vtkActor* actor)
{
    vtkMapper* mapper = actor->GetMapper();
    treeManager_->AddDataObject(mapper->GetInput());
}

void IceTCompositer::redistributeActorPolyData(rtppActor& actor)
{
    vtkSmartPointer<vtkPolyData> redistributed = redistributePolyData(actor.getProcessedInputData());
    actor.setProcessedInputData(redistributed);
}

void IceTCompositer::redistributeWidgetActorPolyData(vtkActor* actor)
{
    vtkMapper* mapper = actor->GetMapper();
    vtkPolyData* data = vtkPolyData::SafeDownCast(mapper->GetInput());
    if (data && data->GetNumberOfCells() > 0)
    {
        mapper->SetInputDataObject(redistributePolyData(data));
    }
}

vtkSmartPointer<vtkPolyData> IceTCompositer::redistributePolyData(vtkPolyData* data)
{
    if (!redistributor_->GetPKdTree())
    {
        redistributor_->SetPKdTree(getKdTree());
    }
    redistributor_->SetInputData(data);
    redistributor_->Update();
    vtkNew<vtkPolyDataNormals> normals;
    normals->SetInputData(redistributor_->GetPolyDataOutput());
    normals->Update();
    return normals->GetOutput();
}

void IceTCompositer::renderAndWriteImage
(
        vtkRenderWindow* window,
        const fileName& pathToOutputFileWithoutExtension,
        double currentTime
)
{
    vtkNew<vtkPartitionOrderingInterface> partitionOrderingInterface;
    partitionOrderingInterface->SetImplementation(getKdTree());

    synchronizer_->SetPartitionOrdering(partitionOrderingInterface);
    synchronizer_->SetRenderer(renderer_);

    vtkSmartPointer<vtkSynchronizedRenderWindows> windowSync = createWindowSynchronizer(window);
    renderAndWriteImage_(window, pathToOutputFileWithoutExtension, currentTime);
}


vtkSmartPointer<vtkPKdTree> IceTCompositer::getKdTree()
{
    if (!pKdTree_.Get())
    {
        treeManager_->GenerateKdTree();
        pKdTree_ = treeManager_->GetKdTree();
    }
    return pKdTree_;
}

} // End namespace

// ************************************************************************* //
#endif