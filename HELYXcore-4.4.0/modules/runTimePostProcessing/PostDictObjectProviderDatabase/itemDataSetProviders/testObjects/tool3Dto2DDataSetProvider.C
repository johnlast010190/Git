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

#include <utility>

#include "tool3Dto2DDataSetProvider.H"
#include "vtkMultiProcessController.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkAppendPolyData.h"

#include "engys3Dto2DTool.h"


namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Tool3Dto2DDataSetProvider::Tool3Dto2DDataSetProvider(const std::string& name, Tool3Dto2DObjectData dictData)
: ItemDataSetProvider(name), data_(std::move(dictData))
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Tool3Dto2DDataSetProvider::update(scalar vtkNotUsed(currentTime))
{
    vtkNew<engys3Dto2DTool> tool;
    tool->SetPlaneNormal(data_.normal[0], data_.normal[1], data_.normal[2]);
    tool->SetPlaneOrigin(data_.origin[0], data_.origin[1], data_.origin[2]);
    tool->SetController(vtkMultiProcessController::GetGlobalController());
    for (ItemDataSetProvider* provider : sources_)
    {
        tool->AddInputDataObject(provider->getDataSetOutput());
    }
    tool->Update();
    vtkMultiBlockDataSet* allRibbonSurfaces = tool->GetOutput();

    if (allRibbonSurfaces->GetNumberOfBlocks() == 0)
    {
        output_ = vtkSmartPointer<vtkPolyData>::New();
    }
    else if (allRibbonSurfaces->GetNumberOfBlocks() == 1)
    {
        output_ = allRibbonSurfaces->GetBlock(0);
    }
    else
    {
        vtkNew<vtkAppendPolyData> append;
        for (unsigned int block = 0; block < allRibbonSurfaces->GetNumberOfBlocks(); block++)
        {
            append->AddInputData(vtkPolyData::SafeDownCast(allRibbonSurfaces->GetBlock(block)));
        }
        append->Update();
        output_ = append->GetOutput();
    }
}

} // End namespace


// ************************************************************************* //
