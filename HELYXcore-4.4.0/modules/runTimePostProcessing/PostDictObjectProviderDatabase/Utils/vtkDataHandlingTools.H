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
    (c) 2020-2023 Engys Ltd.

Description
    A collection of static helper functions for dealing with VTK data.

SourceFiles
    vtkDataHandlingTools.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "polySurface/fields/polySurfaceFields.H"

#include "vtkSmartPointer.h"

class vtkPolyData;
class vtkDataSet;
class vtkDataObject;
class vtkMultiProcessController;
class vtkMultiBlockDataSet;

namespace Foam
{
class fvMesh;
namespace vtk::Tools
{

vtkSmartPointer<vtkPolyData> allGatherPolyDataFromProcesses
(
    vtkMultiProcessController* Controller,
    vtkPolyData* dataset
);

vtkSmartPointer<vtkPolyData> gatherPolyDataOnProcess
(
    vtkMultiProcessController* Controller,
    vtkPolyData* dataset,
    int processId
);


// Interpolate cell to point data, returning a data set containing both cell and
// point data.
vtkSmartPointer<vtkDataSet> convertCellToPoint(vtkDataSet* inputPort);

vtkSmartPointer<vtkPolyData> dataObjectToPolyData(vtkDataObject* inputData);

[[nodiscard]] vtkSmartPointer<vtkDataSet> removeGhostCells(vtkDataSet *data);

[[nodiscard]] vtkSmartPointer<vtkMultiBlockDataSet> readSurfaceFile(const fileName& filePath);

[[nodiscard]] vtkSmartPointer<vtkDataSet> mergeMultiBlockDataSet(vtkMultiBlockDataSet* multiBlockDataSet);

} // End namespace vtk
} // End namespace Foam
