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

\*---------------------------------------------------------------------------*/

#include "streamlineSourceProviders.H"

#include "Utils/vtkDataHandlingTools.H"

#include <utility>
#include "uniformDistribution/gatheredDomainUniformInput.H"

#include "vtkPolyData.h"
#include "vtkTriangleFilter.h"
#include "vtkMassProperties.h"
#include "vtkPolyDataPointSampler.h"
#include "vtkMultiProcessController.h"

namespace Foam::functionObjects::runTimeVis
{

PatchSourceProvider::PatchSourceProvider(PatchSourceData  data)
:
        data_(std::move(data))
{
}

vtkSmartPointer<vtkPolyData> PatchSourceProvider::getSeedPoints()
{
    vtkMultiProcessController* controller = vtkMultiProcessController::GetGlobalController();

    vtkSmartPointer<vtkPolyData> fullPatchDataSet;

    if (controller)
    {
        fullPatchDataSet =
            vtk::Tools::allGatherPolyDataFromProcesses
            (
                controller,
                patchSource_
            );
    }
    else
    {
        fullPatchDataSet = patchSource_;
    }

    vtkIdType patchSourceNumberOfPoints = fullPatchDataSet->GetNumberOfPoints();

    if (data_.numberOfPoints < patchSourceNumberOfPoints) {
        int seedNumber = 10000;
        GatheredDomainUniformInput uniformInput
        (
                fullPatchDataSet,
                data_.numberOfPoints,
                seedNumber
        );
        return uniformInput.compute();
    } else {
        vtkNew<vtkTriangleFilter> triangleFilter;
        triangleFilter->SetInputData(fullPatchDataSet);
        vtkNew<vtkMassProperties> massFilter;
        massFilter->SetInputConnection(triangleFilter->GetOutputPort());
        massFilter->Update();

        auto area = static_cast<scalar>(massFilter->GetSurfaceArea());
        scalar distance = std::sqrt(area / static_cast<scalar>(data_.numberOfPoints));

        vtkNew<vtkPolyDataPointSampler> sample;
        sample->SetInputData(fullPatchDataSet);
        sample->SetDistance(distance);
        sample->Update();

        return sample->GetOutput();
    }
}

bool PatchSourceProvider::acceptsPatchSource() const
{
    return true;
}

void PatchSourceProvider::setPatchSource(vtkDataObject* patch)
{
    patchSource_ = vtkPolyData::SafeDownCast(patch);

    if (patchSource_ == nullptr) {
        FatalError << "Could not process patch streamline source " << data_.patchId.name << exit(FatalError);
    }
}

} // End namespace

// ************************************************************************* //
