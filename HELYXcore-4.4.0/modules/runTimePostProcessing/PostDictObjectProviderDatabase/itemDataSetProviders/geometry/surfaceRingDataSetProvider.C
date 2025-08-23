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

#include "surfaceRingDataSetProvider.H"

#include "vtkLineSource.h"
#include "vtkTubeFilter.h"
#include "vtkAppendPolyData.h"
#include "vtkTriangleFilter.h"
#include "vtkPolyDataNormals.h"

namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void SurfaceRingDataSetProvider::initialiseSurface()
{
    if (!ParallelUtils::isMaster())
    {
        output_ = vtkSmartPointer<vtkPolyData>::New();
        return;
    }

    vtkNew<vtkLineSource> lineSource;
    lineSource->SetPoint1(data_.point1.x(), data_.point1.y(), data_.point1.z());
    lineSource->SetPoint2(data_.point2.x(), data_.point2.y(), data_.point2.z());

    // Create a tube (cylinder) around the line
    vtkNew<vtkTubeFilter> internalTubeFilter;
    internalTubeFilter->SetInputConnection(lineSource->GetOutputPort());
    internalTubeFilter->SetCapping(0);
    internalTubeFilter->SetRadius(data_.innerRadius);
    internalTubeFilter->SetNumberOfSides(50);
    // outputData_ = internalTubeFilter->GetOutput();

    vtkNew<vtkTubeFilter> externalTubeFilter;
    externalTubeFilter->SetInputConnection(lineSource->GetOutputPort());
    externalTubeFilter->SetCapping(0);
    externalTubeFilter->SetRadius(data_.outerRadius);
    externalTubeFilter->SetNumberOfSides(50);

    vtkNew<vtkAppendPolyData> append;
    append->AddInputConnection(internalTubeFilter->GetOutputPort());
    append->AddInputConnection(externalTubeFilter->GetOutputPort());
    append->Update();

    vtkSmartPointer<vtkPolyData> outputMesh = append->GetOutput();
    // Now cap the ring

    // Don't triangulate if point1 and point2 are equal, otherwise you get
    // an exception in dataset.getBounds()
    if(data_.point1 != data_.point2)
    {
        auto filter = vtkSmartPointer<vtkTriangleFilter>::New();
        filter->SetInputData(outputMesh);
        filter->Update();
        outputMesh = filter->GetOutput();
    }

    vtkSmartPointer<vtkCellArray> outputTriangles = outputMesh->GetPolys();

    vtkIdType length = internalTubeFilter->GetOutput()->GetNumberOfPoints();
    for (int ptId = 0; ptId < 50; ptId++)
    {
        // Triangle one extremity
        auto triangle = vtkSmartPointer<vtkIdList>::New();
        triangle->InsertNextId(ptId);
        triangle->InsertNextId(ptId + length);
        triangle->InsertNextId((ptId + 1) % 50 + length);
        outputTriangles->InsertNextCell(triangle);

        triangle = vtkSmartPointer<vtkIdList>::New();
        triangle->InsertNextId(ptId);
        triangle->InsertNextId((ptId + 1) % 50 + length);
        triangle->InsertNextId((ptId + 1) % 50);
        outputTriangles->InsertNextCell(triangle);

        // Triangle the other extremity
        vtkIdType offset = length - 50;
        triangle = vtkSmartPointer<vtkIdList>::New();
        triangle->InsertNextId(ptId + offset);
        triangle->InsertNextId(ptId + +offset + length);
        triangle->InsertNextId((ptId + 1) % 50 + offset + length);
        outputTriangles->InsertNextCell(triangle);

        triangle = vtkSmartPointer<vtkIdList>::New();
        triangle->InsertNextId((ptId + 1) % 50 + length + offset);
        triangle->InsertNextId((ptId + 1) % 50 + offset);
        triangle->InsertNextId(ptId + offset);
        outputTriangles->InsertNextCell(triangle);
    }

    vtkNew<vtkPolyDataNormals> normals;
    normals->SetInputData(outputMesh);
    normals->Update();
    output_ = normals->GetOutput();
}

} // End namespace

// ************************************************************************* //
