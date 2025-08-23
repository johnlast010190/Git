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
    (c) 2024-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include <utility>

#include "streamlinesDataSetProvider.H"
#include "Utils/vtkDataHandlingTools.H"
#include "Utils/ParallelUtils.H"

#include "vtkMinimalStandardRandomSequence.h"
#include "vtkStreamTracer.h"
#include "vtkPStreamTracer.h"
#include "vtkTubeFilter.h"
#include "vtkPointData.h"
#include "vtkCellData.h"


namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

StreamlinesDataSetProvider::StreamlinesDataSetProvider
(
        const std::string& name,
        const StreamlinesObjectData& dictData,
        std::shared_ptr<StreamlineSourceProvider> sourceProvider
)
: ItemDataSetProvider(name),
  sourceProvider_(std::move(sourceProvider))
{
    streamTracer_ = vtkSmartPointer<vtkStreamTracer>::New();

    if (ParallelUtils::isRunningInParallel())
    {
#if !defined( WIN32 ) && !defined ( WIN64 )
        // WARNING:  PStreamTracer must have the same seed data on all processes
        streamTracer_ = vtkSmartPointer<vtkPStreamTracer>::New();
#endif
    }
    streamTracer_->SetInputArrayToProcess
    (
        0,  // index: vectors(0)
        0,  // port
        0,  // connection
        vtkDataObject::FIELD_ASSOCIATION_POINTS,
        dictData.vectorField.c_str()
    );

    if (!sourceProvider_->acceptsPatchSource())
    {
        streamTracer_->SetSourceData(sourceProvider_->getSeedPoints());
    }
    streamTracer_->SetInterpolatorTypeToDataSetPointLocator();
    streamTracer_->SetIntegrationDirectionToBoth();
    streamTracer_->SetIntegratorTypeToRungeKutta45();
    streamTracer_->SetInitialIntegrationStep(0.2);
    streamTracer_->SetMinimumIntegrationStep(0.01);
    streamTracer_->SetMaximumIntegrationStep(0.5);
    streamTracer_->SetComputeVorticity(false);
    streamTracer_->SetTerminalSpeed(1e-12);
    streamTracer_->SetMaximumError(1e-6);
    streamTracer_->SetMaximumNumberOfSteps(dictData.maxSteps);
    streamTracer_->SetMaximumPropagation(dictData.maxLength);

    tubeFilter_ = vtkSmartPointer<vtkTubeFilter>::New();
    tubeFilter_->SetInputConnection(streamTracer_->GetOutputPort());
    tubeFilter_->SetInputArrayToProcess
    (
        0,  // index: vectors(0)
        0,  // port
        0,  // connection
        vtkDataObject::FIELD_ASSOCIATION_POINTS,
        dictData.vectorField.c_str()
    );
    tubeFilter_->SetRadius(dictData.radius);
    tubeFilter_->SetVaryRadiusToVaryRadiusOff();
    tubeFilter_->SetNumberOfSides(12);
    tubeFilter_->CappingOn();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void StreamlinesDataSetProvider::update(scalar vtkNotUsed(currentTime))
{
    streamTracer_->SetInputData(sources_.at(0)->getDataSetOutput());
    if (sourceProvider_->acceptsPatchSource())
    {
        sourceProvider_->setPatchSource(sources_.at(1)->getDataSetOutput());
        vtkSmartPointer<vtkPolyData> seedPoints = sourceProvider_->getSeedPoints();
        streamTracer_->SetSourceData(seedPoints);
    }
    tubeFilter_->Update();
    vtkSmartPointer<vtkDataSet> output = tubeFilter_->GetOutput();
    if (output->GetNumberOfPoints() > 0)
    {
        output->GetCellData()->RemoveArray("ReasonForTermination");
        output->GetCellData()->RemoveArray("SeedIds");
        output->GetPointData()->RemoveArray("IntegrationTime");
        output->GetPointData()->RemoveArray("TubeNormals");
    }
    output_ = output;
}

} // End namespace


// ************************************************************************* //
