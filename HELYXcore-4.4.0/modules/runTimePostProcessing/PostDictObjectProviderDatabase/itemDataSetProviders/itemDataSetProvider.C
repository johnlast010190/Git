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
    (c) 2022-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include <utility>

#include "itemDataSetProvider.H"
#include "vtkMultiBlockDataSet.h"
#include "engysOutlineFilter.h"
#include "vtkMultiProcessController.h"

#include "Utils/basicProfiler.H"

#undef Log
#include "vtkLogger.h"


namespace Foam::functionObjects::runTimeVis
{

ItemDataSetProvider::~ItemDataSetProvider() = default;

ItemDataSetProvider::ItemDataSetProvider(std::string  name)
:
    itemRequirements_(),
    name_(std::move(name))
{
    latestProcessedTimeIndex_ = -1e9;
    nonExistentRequiredFields_.clear();
}

void ItemDataSetProvider::addItemRequirements(const ItemRequirements& itemRequirements)
{
    itemRequirements_.mergeFromDownstream(itemRequirements);
    addExtraRequirementsDueToUpstreamRequirements(itemRequirements_, itemRequirements);
}

void ItemDataSetProvider::addToNonExistingFields(const foamField& field)
{
    if (!nonExistentRequiredFields_.containsFoamField(field))
    {
        WarningInFunction << "Did not find " << field
                          << " as a vector or scalar field for "
                          << name_ << endl;
        nonExistentRequiredFields_.addField(field);
    }
}

void ItemDataSetProvider::addItemRequirementsToSources() {
    for (ItemDataSetProvider *source : sources_) {
        source->addItemRequirementsForThisAndSources(itemRequirements_);
    }
}

void ItemDataSetProvider::addItemRequirementsForThisAndSources(const ItemRequirements& itemRequirements)
{
    addItemRequirements(itemRequirements);
    addItemRequirementsToSources();
}

void ItemDataSetProvider::setInstanceServicesForThisAndSources(InstanceServices *instanceServices)
{
    this->instanceServices_ = instanceServices;
    for (ItemDataSetProvider* source : sources_)
    {
        source->setInstanceServicesForThisAndSources(instanceServices);
    }
}

void ItemDataSetProvider::updateIfNecessary(Foam::label currentTimeIndex, Foam::scalar currentTime)
{
    if (needsToUpdate(currentTimeIndex))
    {
        latestProcessedTimeIndex_ = currentTimeIndex;
        for (ItemDataSetProvider* source : getUsedSources())
        {
            source->updateIfNecessary(currentTimeIndex, currentTime);
        }
        vtkLogF(INFO, "Updating %s", name_.c_str());
#if BASIC_PROFILE == 1
        basicProfiler::startGroup("updating " + name_);
#endif
        update(currentTime);
#if BASIC_PROFILE == 1
        basicProfiler::addMeasurePoint("finished updating " + name_);
        basicProfiler::endGroup();
#endif
    }
}

vtkSmartPointer<vtkDataSet> ItemDataSetProvider::getDataSetOutput(unsigned int fallbackBlock)
{
    vtkDataSet* dataset = vtkDataSet::SafeDownCast(output_);
    if (dataset) return dataset;

    vtkMultiBlockDataSet *multiBlockInput = vtkMultiBlockDataSet::SafeDownCast(output_);
    if (multiBlockInput && multiBlockInput->GetNumberOfBlocks() > fallbackBlock)
    {
        return vtkDataSet::SafeDownCast(multiBlockInput->GetBlock(fallbackBlock));
    }
    else
    {
        return vtkSmartPointer<vtkPolyData>::New();
    }
}

vtkSmartPointer<vtkPolyData> ItemDataSetProvider::getDataSetOutline()
{
    vtkSmartPointer<vtkDataSet> dataSet = getDataSetOutput();

    vtkNew<vtkFieldData> inputFieldData;
    inputFieldData->DeepCopy(dataSet->GetFieldData());

    vtkNew<engysOutlineFilter> outlineFilter;
    outlineFilter->AllReduceOff();
    outlineFilter->SetController(vtkMultiProcessController::GetGlobalController());

    outlineFilter->SetInputData(dataSet);
    outlineFilter->Update();
    vtkSmartPointer<vtkPolyData> outline = outlineFilter->GetOutput();
    outline->SetFieldData(inputFieldData);
    return outline;
}

} // End namespace Foam


// ************************************************************************* //
