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

#include "cellZoneDataSetProvider.H"

#include "vtk/adaptor/foamVtkVtuAdaptor.H"
#include "fvMeshSubsetProxy/fvMeshSubsetProxy.H"

#include "engysExtractCells.h"

namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void CellZoneDataSetProvider::update(scalar currentTime)
{

    const zone* cellZoneObj = findCellZoneObject();

    if (!cellZoneObj) {
        vtkNew<vtkUnstructuredGrid> empty;
        output_ = setGhostCells(empty);
        return;
    }

    output_ = setGhostCells(extractCellZoneObject(cellZoneObj));
}

const zone* CellZoneDataSetProvider::findCellZoneObject()
{
    const cellZoneMesh *cellZoneMeshObj;
    try
    {
        cellZoneMeshObj = &meshAndFields_.getMesh().cellZones();
    }
    catch (const std::exception &ex)
    {
        FatalErrorInFunction << "Error reading the cell zones: "
                             << ex.what() << abort(FatalError);
        return nullptr;
    }

    const zone *cellZoneObj;
    try
    {
        cellZoneObj = &(*cellZoneMeshObj)[getName()];
    }
    catch (const std::exception &ex)
    {
        FatalErrorInFunction << "Error reading cell zone " << getName() << ": "
                             << ex.what() << abort(FatalError);
        return nullptr;
    }
    return cellZoneObj;
}

vtkSmartPointer<vtkDataSet> CellZoneDataSetProvider::extractCellZoneObject(const zone* cellZone)
{
    vtkNew<vtkIdTypeArray> cellIdList;
    cellIdList->SetNumberOfComponents(1);
    label nCells = cellZone->size();
    cellIdList->SetNumberOfValues(nCells);
    for (int i = 0; i < nCells; i++)
    {
        cellIdList->SetValue(i, cellZone->operator[](i));
    }
    vtkNew<engysExtractCells> cellExtractor;
    cellExtractor->SetInputData(sources_.at(0)->getDataSetOutput());
    cellExtractor->SetCellIdList(cellIdList);
    cellExtractor->Update();

    vtkSmartPointer<vtkDataSet> convertedMesh = cellExtractor->GetOutput();

    convertedMesh->GetCellData()->RemoveArray("vtkOriginalCellIds");
    if (convertedMesh->GetNumberOfPoints() > 0)
    {
        convertedMesh->GetPointData()->GetArray("vtkOriginalPointIds")->SetName("engysGhostMeshPointsArray");
        convertedMesh->GetPointData()->Modified();
    }
    return convertedMesh;
}

} // End namespace Foam


// ************************************************************************* //
