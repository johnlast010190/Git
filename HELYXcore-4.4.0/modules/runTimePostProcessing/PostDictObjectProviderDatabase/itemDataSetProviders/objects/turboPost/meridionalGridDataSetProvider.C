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
    (c) 2023-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include <vtkAppendPolyData.h>
#include "meridionalGridDataSetProvider.H"

#include "dataStructs/objects/turboPost/meridionalGridObjectData.H"

#include "engysTurboPostBase.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkMultiProcessController.h"
#include "engysMeridionalPlot.h"

namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

MeridionalGridDataSetProvider::MeridionalGridDataSetProvider(
    const std::string &name,
    const MeridionalGridObjectData &dictData
)
    : ItemDataSetProvider(name)
{
    this->inletPatches = dictData.inletPatches.size();
    this->outletPatches = dictData.outletPatches.size();
    this->hubPatches = dictData.hubPatches.size();
    this->shroudPatches = dictData.shroudPatches.size();

    this->turboBackgroundMesh_ = vtkSmartPointer<engysTurboPostBase>::New();
    this->turboBackgroundMesh_->SetSpanwiseElements(dictData.spanElements);
    this->turboBackgroundMesh_->SetStreamwiseElements(dictData.streamElements);
    this->turboBackgroundMesh_->SetOrigin(dictData.origin[0], dictData.origin[1], dictData.origin[2]);
    this->turboBackgroundMesh_->SetAxis(dictData.axis[0], dictData.axis[1], dictData.axis[2]);
    this->turboBackgroundMesh_->SetXDirection(dictData.xDirection[0], dictData.xDirection[1], dictData.xDirection[2]);

    vtkMultiProcessController* controller = vtkMultiProcessController::GetGlobalController();
    this->turboBackgroundMesh_->SetController(controller);

    this->meridionalPlot_ = vtkSmartPointer<engysMeridionalPlot>::New();
    this->meridionalPlot_->SetController(controller);
    this->meridionalPlot_->SetInputConnection(turboBackgroundMesh_->GetOutputPort());
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void MeridionalGridDataSetProvider::update(scalar currentTime)
{
    int sourceIndex = 0;
    vtkDataSet *domainDataSet = sources_.at(sourceIndex++)->getDataSetOutput();
    turboBackgroundMesh_->SetInputData(domainDataSet);

    turboBackgroundMesh_->RemoveAllInletInputs();
    for (int i = 0; i < inletPatches; i++)
    {
        vtkPolyData *inlet = vtkPolyData::SafeDownCast(sources_.at(sourceIndex++)->getDataSetOutput());
        turboBackgroundMesh_->AddInletInput(inlet);
    }

    turboBackgroundMesh_->RemoveAllOutletInputs();
    for (int i = 0; i < outletPatches; i++)
    {
        vtkPolyData *inlet = vtkPolyData::SafeDownCast(sources_.at(sourceIndex++)->getDataSetOutput());
        turboBackgroundMesh_->AddOutletInput(inlet);
    }

    turboBackgroundMesh_->RemoveAllHubInputs();
    for (int i = 0; i < hubPatches; i++)
    {
        vtkPolyData *inlet = vtkPolyData::SafeDownCast(sources_.at(sourceIndex++)->getDataSetOutput());
        turboBackgroundMesh_->AddHubInput(inlet);
    }

    turboBackgroundMesh_->RemoveAllShroudInputs();
    for (int i = 0; i < shroudPatches; i++)
    {
        vtkPolyData *inlet = vtkPolyData::SafeDownCast(sources_.at(sourceIndex++)->getDataSetOutput());
        turboBackgroundMesh_->AddShroudInput(inlet);
    }

    meridionalPlot_->Update();
    output_ = meridionalPlot_->GetOutput();
}

} // End namespace Foam


// ************************************************************************* //
