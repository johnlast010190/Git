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

#include "surfaceSplitterDataSetProvider.H"

#include "engysSplitByFeatureAngle.h"
#include "vtkMultiProcessController.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkFloatArray.h"
#include "vtkCellData.h"
#include "Utils/vtkDataHandlingTools.H"


namespace Foam::functionObjects::runTimeVis
{

SurfaceSplitterDataSetProvider::SurfaceSplitterDataSetProvider
(
        const std::string& name,
        const SurfaceSplitterObjectData& dictData
)
:
        ItemDataSetProvider(name),
        data_(dictData)
{
}

void SurfaceSplitterDataSetProvider::update(scalar vtkNotUsed(currentTime))
{
    vtkNew<engysSplitByFeatureAngle> splitter;
    splitter->SetController(vtkMultiProcessController::GetGlobalController());

    splitter->SetInputDataObject(this->sources_[0]->getDataObjectOutput());
    splitter->SetAngleDegrees(data_.angleDegrees);
    splitter->Update();

    vtkSmartPointer<vtkMultiBlockDataSet> mb = splitter->GetOutput();
    for (unsigned int i = 0; i < mb->GetNumberOfBlocks(); i++)
    {
        vtkDataSet* block = vtkDataSet::SafeDownCast(mb->GetBlock(i));
        vtkNew<vtkFloatArray> idCellArray;
        idCellArray->SetNumberOfComponents(1);
        idCellArray->SetNumberOfValues(block->GetNumberOfCells());
        idCellArray->SetName("p");
        idCellArray->Fill(i);
        block->GetCellData()->AddArray(idCellArray);
    }

    output_ = vtk::Tools::mergeMultiBlockDataSet(mb);
}

} // End namespace


// ************************************************************************* //
