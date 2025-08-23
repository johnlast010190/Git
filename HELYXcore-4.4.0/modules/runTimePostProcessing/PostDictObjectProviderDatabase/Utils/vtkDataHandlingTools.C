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
    (c) 2020-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "vtkDataHandlingTools.H"

#include "Utils/ParallelUtils.H"

#include "vtkMultiProcessController.h"
#include "vtkCompositeDataGeometryFilter.h"
#include "vtkMultiPieceDataSet.h"
#include "vtkDataSetSurfaceFilter.h"
#include "engysParallelUtils.h"
#include "engysRemoveGhostCells.h"
#include "vtkMultiBlockDataSet.h"
#include "engysSTLReader.h"
#include "vtkBTSReader.h"
#include "engysOBJReader.h"
#include "engysNASReader.h"
#include "engysCADReader.h"
#include "vtkAppendPolyData.h"
#include "vtkAppendFilter.h"
#include "vtkUnstructuredGrid.h"
#include "types/surfaceFormat.H"

namespace Foam::vtk::Tools
{

vtkSmartPointer<vtkPolyData> allGatherPolyDataFromProcesses
(
    vtkMultiProcessController* Controller,
    vtkPolyData* dataset
)
{
    // We're going to collect on proc 0.  It doesn't really matter whether
    // that's the root proc or not, but we need to be consistent throughout this
    // method.
    vtkSmartPointer<vtkPolyData> output = gatherPolyDataOnProcess(Controller, dataset, 0);
    Controller->Broadcast(output, 0);
    return output;
}

vtkSmartPointer<vtkPolyData> gatherPolyDataOnProcess
(
    vtkMultiProcessController* controller,
    vtkPolyData* dataset,
    int processId
)
{
    if (!controller || controller->GetNumberOfProcesses() <= 1)
    {
        return dataset;
    }
    vtkNew<engysParallelUtils>utils;
    utils->SetController(controller);
    return utils->gatherPolyDataOnProcess(dataset, processId);
}


vtkSmartPointer<vtkPolyData> dataObjectToPolyData(vtkDataObject* inputData)
{
    if (inputData->IsA("vtkMultiPieceDataSet"))
    {
        vtkMultiPieceDataSet* mpds = vtkMultiPieceDataSet::SafeDownCast(inputData);
        inputData = mpds->GetPiece(functionObjects::runTimeVis::ParallelUtils::localProcessId());
    }

    if (inputData->IsA("vtkPolyData"))
    {
        return vtkPolyData::SafeDownCast(inputData);
    }
    else
    {
        vtkNew<vtkDataSetSurfaceFilter> filterClipper;
        filterClipper->SetInputData(inputData);
        filterClipper->Update();
        return filterClipper->GetOutput();
    }
}

vtkSmartPointer<vtkDataSet> removeGhostCells(vtkDataSet *data)
{
    vtkNew<engysRemoveGhostCells> ghostCellsRemover;
    ghostCellsRemover->SetInputData(data);
    ghostCellsRemover->Update();
    return ghostCellsRemover->GetOutput();
}

vtkSmartPointer<vtkMultiBlockDataSet> readSurfaceFile(const fileName& filePath)
{
    functionObjects::runTimeVis::SurfaceFormat surfaceFormat =
        functionObjects::runTimeVis::SurfaceFormat::fromNameWithExtension(filePath.name());

    switch(surfaceFormat.getValue())
    {

        case functionObjects::runTimeVis::SurfaceFormat::STL_ASCII:
        case functionObjects::runTimeVis::SurfaceFormat::STL_BINARY:
        case functionObjects::runTimeVis::SurfaceFormat::STL_BINARY_ENGYS:
        case functionObjects::runTimeVis::SurfaceFormat::STL_GZ_ASCII:
        {
            vtkNew<engysSTLReader> reader;
            reader->SetFileName(filePath.c_str());
            reader->Update();
            return reader->GetOutput();
        }

        case functionObjects::runTimeVis::SurfaceFormat::OBJ:
        {
            vtkNew<engysOBJReader> reader;
            reader->SetFileName(filePath.c_str());
            reader->Update();
            return reader->GetOutput();
        }

        case functionObjects::runTimeVis::SurfaceFormat::BTS:
        {
            vtkNew<vtkBTSReader> reader;
            reader->SetFileName(filePath.c_str());
            reader->Update();
            return reader->GetOutput();
        }

        case functionObjects::runTimeVis::SurfaceFormat::NAS:
        case functionObjects::runTimeVis::SurfaceFormat::NAS_GZ:
        {
            vtkNew<engysNASReader> reader;
            reader->SetFileName(filePath.c_str());
            reader->Update();
            return reader->GetOutput();
        }
        case functionObjects::runTimeVis::SurfaceFormat::STP:
        {
            vtkNew<engysCADReader> reader;
            reader->SetFileName(filePath.c_str());
            reader->SetFileFormat(vtkOCCTReader::STEP);
            reader->Update();
            return reader->GetOutput();
        }

        default:
            Warning
                << "The file has an unknown extension: " << filePath.ext()
                << endl;
            return nullptr;
    }
}

static bool isAllPolyData(vtkMultiBlockDataSet* multiBlockDataSet)
{
    for (unsigned int blockId = 0; blockId < multiBlockDataSet->GetNumberOfBlocks(); blockId++)
    {
        vtkDataObject* block = multiBlockDataSet->GetBlock(blockId);
        if (vtkPolyData::SafeDownCast(block)) continue;
        if (vtkMultiBlockDataSet* subMB = vtkMultiBlockDataSet::SafeDownCast(block))
        {
            if (!isAllPolyData(subMB))
            {
                return false;
            }
        }
        else
        {
            return false;
        }
    }
    return true;
}

static void addDataSets(vtkAlgorithm* algorithm, vtkMultiBlockDataSet* multiBlockDataSet)
{
    for (unsigned int blockId = 0; blockId < multiBlockDataSet->GetNumberOfBlocks(); blockId++)
    {
        vtkDataObject* block = multiBlockDataSet->GetBlock(blockId);
        if (vtkDataSet::SafeDownCast(block))
        {
            algorithm->AddInputDataObject(block);
        }
        else if (vtkMultiBlockDataSet* subMB = vtkMultiBlockDataSet::SafeDownCast(block))
        {
            addDataSets(algorithm, subMB);
        }
    }
}

vtkSmartPointer<vtkDataSet> mergeMultiBlockDataSet(vtkMultiBlockDataSet* multiBlockDataSet)
{
    if (multiBlockDataSet->GetNumberOfBlocks() == 0)
    {
        return vtkSmartPointer<vtkPolyData>::New();
    }
    else if (multiBlockDataSet->GetNumberOfBlocks() == 1)
    {
        if (vtkDataSet* dataSet = vtkDataSet::SafeDownCast(multiBlockDataSet->GetBlock(0)))
        {
            return dataSet;
        }
        if (vtkMultiBlockDataSet* subMB = vtkMultiBlockDataSet::SafeDownCast(multiBlockDataSet->GetBlock(0)))
        {
            return mergeMultiBlockDataSet(subMB);
        }
        return nullptr;
    }
    else if (isAllPolyData(multiBlockDataSet))
    {
        vtkNew<vtkAppendPolyData> append;
        addDataSets(append, multiBlockDataSet);
        append->Update();
        return append->GetOutput();
    }
    else
    {
        vtkNew<vtkAppendFilter> append;
        addDataSets(append, multiBlockDataSet);
        append->Update();
        return append->GetOutput();
    }
}

} // End namespace Foam
