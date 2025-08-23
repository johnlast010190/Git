/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.3.0
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
    (c) 2025 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#include "rtppVtkVtuAdaptor.H"
#include "engysVectorCoreDataArray.H"

namespace Foam::functionObjects::runTimeVis
{

vtkSmartPointer<vtkPoints> rtppVtuAdaptor::points(const Foam::fvMesh &mesh) const
{
    // Convert OpenFOAM mesh vertices to VTK

    // Normal points
    const pointField &pts = mesh.points();

    vtkNew<engysVectorCoreDataArray> wrappedPoints;
    wrappedPoints->SetCorePointer(&pts);

    // Additional cell centres
    const labelUList &addPoints = this->adaptor.additionalIds();

    wrappedPoints->SetNumberOfTuples(pts.size() + addPoints.size());

    vtkIdType pointId = pts.size();
    // Cell centres
    for (const label meshCelli: addPoints)
    {
        wrappedPoints->SetTypedTuple(pointId++, mesh.cellCentres()[meshCelli].v_);
    }

    vtkNew<vtkPoints> vtkpoints;
    vtkpoints->SetData(wrappedPoints);

    return vtkpoints;
}


vtkSmartPointer<vtkUnstructuredGrid> rtppVtuAdaptor::internal
    (
        const fvMesh &mesh,
        const bool decompPoly
    ) const
{
    vtk::vtuSizing sizing(mesh, decompPoly);

    auto cellTypes = vtkSmartPointer<vtkUnsignedCharArray>::New();

    auto cells = vtkSmartPointer<vtkCellArray>::New();
    auto faces = vtkSmartPointer<vtkIdTypeArray>::New();

    auto cellLocations = vtkSmartPointer<vtkIdTypeArray>::New();
    auto faceLocations = vtkSmartPointer<vtkIdTypeArray>::New();

    UList<uint8_t> cellTypesUL =
        vtk::Tools::asUList(cellTypes, sizing.nFieldCells());

    List<vtkIdType> cellsL(sizing.sizeInternal(vtk::vtuSizing::slotType::CELLS));
    cells->AllocateEstimate(sizing.nFieldCells(), 10);

    UList<vtkIdType> cellLocationsUL =
        vtk::Tools::asUList
            (
                cellLocations,
                sizing.sizeInternal(vtk::vtuSizing::slotType::CELLS_OFFSETS)
            );

    UList<vtkIdType> facesUL =
        vtk::Tools::asUList
            (
                faces,
                sizing.sizeInternal(vtk::vtuSizing::slotType::FACES)
            );

    UList<vtkIdType> faceLocationsUL =
        vtk::Tools::asUList
            (
                faceLocations,
                sizing.sizeInternal(vtk::vtuSizing::slotType::FACES_OFFSETS)
            );


    sizing.populateInternal
        (
            mesh,
            cellTypesUL,
            cellsL,
            cellLocationsUL,
            facesUL,
            faceLocationsUL,
            static_cast<foamVtkMeshMaps &>(this->adaptor)
        );

    auto iter = cellsL.begin();
    const auto &end = cellsL.end();
    while (iter != end)
    {
        int cellSize = int(*iter);
        iter++;
        cells->InsertNextCell(cellSize);
        for (vtkIdType i = 0; i < cellSize; i++)
        {
            cells->InsertCellPoint(*iter);
            iter++;
        }
    }
    cells->Squeeze();


    vtkNew<vtkUnstructuredGrid> vtkmesh;

    // Convert OpenFOAM mesh vertices to VTK
    // - can only do this *after* populating the decompInfo with cell-ids
    //   for any additional points (ie, mesh cell-centres)
    vtkmesh->SetPoints(this->points(mesh));

    if (!facesUL.empty())
    {
        vtkmesh->SetCells
            (
                cellTypes,
                cells,
                faceLocations,
                faces
            );
    }
    else
    {
        vtkmesh->SetCells
            (
                cellTypes,
                cells,
                nullptr,
                nullptr
            );
    }

    return vtkmesh;
}

}

// ************************************************************************* //
