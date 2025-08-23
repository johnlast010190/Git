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
    (c) 2015-2019 OpenCFD Ltd.
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include <memory>
#include <utility>

#include "instanceServices.H"

#include "engysNewSetGhostCells.h"
#include "postDictObjectProviderDatabaseInstance.H"
#include "vtkMultiProcessController.h"
#include "engysSGCProcessBoundaries.h"

namespace Foam::functionObjects::runTimeVis
{

InstanceServices::InstanceServices(PostDictObjectProviderDatabaseInstance& instance)
: instance_(instance)
{
    latestProcessedTimeIndex_ = -1e9;
    currentTime_ = 0;
}

void InstanceServices::updateIfNecessary(label currentTimeIndex, scalar currentTime)
{
    if (latestProcessedTimeIndex_ != currentTimeIndex)
    {
        latestProcessedTimeIndex_ = currentTimeIndex;
        currentTime_ = currentTime;
        for (auto &region: regions_)
        {
            region.second->updateIfNecessary(currentTimeIndex, currentTime);
        }
    }
}

void InstanceServices::reset()
{
    regions_.clear();
}

InstanceServicesRegion* InstanceServices::getInstanceServicesForRegion(const std::string& region)
{
    auto itr = regions_.find(region);
    if (itr != regions_.end())
    {
        return itr->second.get();
    }
    auto emplaced = regions_.emplace(
        region,
        std::make_unique<InstanceServicesRegion>(
                region,
                instance_
            ));
    InstanceServicesRegion* services = emplaced.first->second.get();
    services->updateIfNecessary(latestProcessedTimeIndex_, currentTime_);
    return services;
}

//void InstanceServices::updateProvidersIfNecessary(const std::vector<ItemDataSetProvider*>& providers)
//{
//    for (ItemDataSetProvider* provider : providers)
//    {
//        provider->updateIfNecessary(latestProcessedTimeIndex_, currentTime_);
//    }
//}


InstanceServicesRegion::InstanceServicesRegion
    (
        std::string region,
        PostDictObjectProviderDatabaseInstance& instance
    ) :
    region_(std::move(region)),
    instance_(instance)
{
    latestProcessedTimeIndex_ = -1e9;
    currentTime_ = 0;
    reset();
}

void InstanceServicesRegion::updateIfNecessary(label currentTimeIndex, scalar currentTime)
{
    if (latestProcessedTimeIndex_ != currentTimeIndex)
    {
        latestProcessedTimeIndex_ = currentTimeIndex;
        currentTime_ = currentTime;
        reset();
    }
}

void InstanceServicesRegion::reset()
{
    this->ghostCells_ = nullptr;
    this->externalMeshBoundaries_.reset(nullptr);
}

engysNewSetGhostCells* InstanceServicesRegion::getInitializedGhostCellsFilter()
{
    if (!this->ghostCells_)
    {
        initializeGhostCellsFilter();
    }
    return this->ghostCells_;
}

void InstanceServicesRegion::initializeGhostCellsFilter()
{
    vtkMultiProcessController* controller = vtkMultiProcessController::GetGlobalController();

    if (controller && controller->GetNumberOfProcesses() > 1)
    {
        const int localProcessId = controller->GetLocalProcessId();

        this->ghostCells_ = vtkSmartPointer<engysNewSetGhostCells>::New();
        ghostCells_->RemoveAllProcessInputs();
        
        auto processProviders = instance_.itemStorage_.getProcessBoundaryProviders(region_);
        
        for (ItemDataSetProvider *boundarySource: processProviders)
        {
            if (boundarySource->isProcessorBoundary())
            {
                boundarySource->updateIfNecessary(latestProcessedTimeIndex_, currentTime_);
                vtkPolyData* boundary = vtkPolyData::SafeDownCast(boundarySource->getDataSetOutput());
                if (!boundary || boundary->GetNumberOfPoints() <= 0 || !boundary->GetFieldData()) continue;

                int boundaryRemoteProcessId = boundarySource->getTargetProcessorBoundary();

                vtkNew<engysSGCProcessBoundaries> extractor;
                extractor->SetLocalRenderProcess(localProcessId);
                extractor->SetRemoteRenderProcess(boundaryRemoteProcessId);
                extractor->AddProcessPatchInput(boundary);
                ghostCells_->SetProcessPatchInput(extractor, boundaryRemoteProcessId);
            }
        }
    }
}

//const std::vector<ItemDataSetProvider*>& InstanceServicesRegion::getExternalBoundaryProviders()
//{
//    if (!this->externalMeshBoundaries_)
//    {
//        this->externalMeshBoundaries_ = std::make_unique<std::vector<ItemDataSetProvider*>>();
//
//        std::vector<ItemDataSetProvider*> allMeshBoundaries = instance_.itemStorage_.getAllMeshBoundaryProviders(region_);
//
//        for (ItemDataSetProvider* meshBoundary : allMeshBoundaries)
//        {
//            if (meshBoundary->isProcessorBoundary())
//            {
//                continue;
//            }
//            else if (!meshBoundary->isNCCPatch() &&
//                     !meshBoundary->isCoupledPatch() &&
//                     !meshBoundary->isEmptyPatch())
//            {
//                externalMeshBoundaries_->push_back(meshBoundary);
//            }
//        }
//    }
//    return *this->externalMeshBoundaries_;
//}

}

// ************************************************************************* //
