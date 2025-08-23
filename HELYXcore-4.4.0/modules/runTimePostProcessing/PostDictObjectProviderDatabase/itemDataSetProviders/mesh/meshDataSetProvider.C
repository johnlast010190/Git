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

#include "meshDataSetProvider.H"

#include "vtkCellDataToPointData.h"
#include "vtkPointDataToCellData.h"
#include "engysNewSetGhostCells.h"
#include "vtkMultiProcessController.h"
#include "engysSGCProcessBoundaries.h"
#include "databaseInstance/instanceServices.H"
#include "meshes/polyMesh/polyPatches/constraint/processor/processorPolyPatch.H"

namespace Foam::functionObjects::runTimeVis
{

MeshDataSetProvider::MeshDataSetProvider(const std::string &name, string region, const MeshAndFields& meshAndFields) :
    ItemDataSetProvider(name),
    region_(std::move(region)),
    meshAndFields_(meshAndFields)
{
}

MeshDataSetProvider::~MeshDataSetProvider() = default;

vtkSmartPointer<vtkDataArray> MeshDataSetProvider::getVolmeshFieldDataArray
    (
        const foamField &fieldName,
        bool& isPointArray
    )
{
    vtkSmartPointer<vtkDataArray> field = meshAndFields_.convertInternalFoamFieldToVtkArray(fieldName,
        isPointArray);
    if (!field)
    {
        addToNonExistingFields(fieldName);
    }
    return field;
}



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
vtkSmartPointer<vtkDataSet> MeshDataSetProvider::setGhostCells(vtkDataSet *mesh, bool isInternalMesh)
{
    bool needGhostCells = getItemRequirements().getNeedsGhostCells();
    bool needPointValues = getRequiredFields().hasPointField();
    if ((!needGhostCells && !needPointValues) || isProcessorBoundary()) return mesh;

    engysNewSetGhostCells* ghostSetter = instanceServices_->getInstanceServicesForRegion(region_)->getInitializedGhostCellsFilter();

    if (ghostSetter)
    {
        ghostSetter->SetInputData(mesh);
        ghostSetter->SetIsInternalMesh(isInternalMesh);
        ghostSetter->SetAddGhostCells(needGhostCells);
        ghostSetter->Update();
        vtkSmartPointer<vtkDataSet> output;
        output.TakeReference(mesh->NewInstance());
        output->ShallowCopy(ghostSetter->GetOutput());
        return output;
    }
    else
    {
        return mesh;
    }
}

void MeshDataSetProvider::addInterpolatedCellAndPointData(vtkDataSet *dataset)
{
    vtkSmartPointer<vtkPointData> interpolatedPointData = interpolateCellToPoint(dataset);
    vtkSmartPointer<vtkCellData> interpolatedCellData = interpolatePointToCell(dataset);
    addInterpolatedData(dataset->GetPointData(), interpolatedPointData);
    addInterpolatedData(dataset->GetCellData(), interpolatedCellData);
}

vtkSmartPointer<vtkCellData> MeshDataSetProvider::interpolatePointToCell(vtkDataSet *dataset)
{
    if (dataset->GetPointData() && dataset->GetPointData()->GetNumberOfArrays() > 0)
    {
        vtkNew<vtkPointDataToCellData> pointToCell;
        pointToCell->PassPointDataOff();
        pointToCell->SetInputData(dataset);
        pointToCell->Update();
        return pointToCell->GetOutput()->GetCellData();
    }
    else
    {
        return vtkSmartPointer<vtkCellData>::New();
    }
}

vtkSmartPointer<vtkPointData> MeshDataSetProvider::interpolateCellToPoint(vtkDataSet *dataset)
{
    if (dataset->GetCellData() && dataset->GetCellData()->GetNumberOfArrays() > 0 && getRequiredFields().hasPointField())
    {
        vtkSmartPointer<vtkDataSet> shallowCopy;
        shallowCopy.TakeReference(dataset->NewInstance());
        shallowCopy->ShallowCopy(dataset);
        std::vector<foamField> cellOnlyFields = getRequiredFields().getCellOnlyFoamFields();
        for (const foamField& cellOnlyField : cellOnlyFields)
        {
            shallowCopy->GetCellData()->RemoveArray(cellOnlyField.c_str());
        }

        vtkNew<vtkCellDataToPointData> cellToPoint;
        cellToPoint->PassCellDataOff();
        cellToPoint->SetInputData(shallowCopy);
        cellToPoint->Update();
        return cellToPoint->GetOutput()->GetPointData();
    }
    else
    {
        return vtkSmartPointer<vtkPointData>::New();
    }
}

void MeshDataSetProvider::addInterpolatedData(vtkDataSetAttributes *baseData, vtkDataSetAttributes *interpolatedData)
{
    for (int i = 0; i < interpolatedData->GetNumberOfArrays(); i++)
    {
        vtkDataArray* array = interpolatedData->GetArray(i);
        baseData->AddArray(array);
    }
}

const FoamFields &MeshDataSetProvider::getRequiredFields()
{
    const FoamFields & requiredFields = getItemRequirements().getRequiredFields();
    if (requiredFields.hasAllFieldsMarker())
    {
        updatedRequiredFields_.clear();
        updatedRequiredFields_.merge(requiredFields);
        updatedRequiredFields_.substituteAllFieldsMarker(meshAndFields_.listAllFields());
        return updatedRequiredFields_;
    }
    return requiredFields;
}

} // End namespace Foam


// ************************************************************************* //
