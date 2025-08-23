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
    (c) 2022-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "surfaceFromFileDataSetProvider.H"

#include "Utils/vtkDataHandlingTools.H"

#include "vtkAppendPolyData.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkInformation.h"
#include "engysMergeCoplanarCells.h"

#undef Log
#include "vtkLogger.h"

namespace Foam::functionObjects::runTimeVis
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void SurfaceFromFileDataSetProvider::initialiseSurface()
{
    if (!ParallelUtils::isMaster())
    {
        output_ = vtkSmartPointer<vtkPolyData>::New();
        return;
    }

    string filePath = data_.caseFolder / "constant" / "triSurface" / data_.fileNameWithExtension;

    // Read in all surfaces from file
    vtkSmartPointer<vtkMultiBlockDataSet> allSurfacesFromFile = readSurfaceFile(filePath);

    output_ = extractDesiredSurface(allSurfacesFromFile);

    if (data_.mergeCoplanar)
    {
        vtkNew<engysMergeCoplanarCells> cellsMerger;
        cellsMerger->SetInputData(output_);
        cellsMerger->Update();
        output_ = cellsMerger->GetOutput();
    }
}

vtkSmartPointer<vtkDataSet> SurfaceFromFileDataSetProvider::extractDesiredSurface(vtkMultiBlockDataSet *allSurfacesFromFile)
{
    if (!allSurfacesFromFile)
    {
        return vtkSmartPointer<vtkPolyData>::New();
    }
    if (data_.readAllSurfaces())
    {
        return vtk::Tools::mergeMultiBlockDataSet(allSurfacesFromFile);
    }
    else
    {
        for (unsigned int block = 0; block < allSurfacesFromFile->GetNumberOfBlocks(); block++)
        {
            std::string blockName = allSurfacesFromFile->GetMetaData(block)->Get(vtkCompositeDataSet::NAME());
            if (blockName == data_.solidSubstring)
            {
                return vtkPolyData::SafeDownCast(allSurfacesFromFile->GetBlock(block));
            }
        }

        Warning << "Solid " << data_.solidSubstring << " was not found in the file!" << endl;

        return vtkSmartPointer<vtkPolyData>::New();
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool isWritable(const fileName& filePath)
{
    autoPtr<ISstream> isPtr(fileHandler().NewIFstream(filePath));
    ISstream& is = isPtr();
    return is.good();
}

vtkSmartPointer<vtkMultiBlockDataSet> SurfaceFromFileDataSetProvider::readSurfaceFile(const fileName& filePath)
{
    if (isWritable(filePath))
    {
        if (data_.readAllSurfaces())
        {
            vtkLogF(INFO, "Reading all surfaces from file %s", filePath.c_str());
        }
        else
        {
            vtkLogF(INFO, "Reading surface %s from file %s", data_.solidSubstring.c_str(), filePath.c_str());
        }
    }
    else
    {
        Warning
            << "The following file could not be read:" << nl
            << "    " << filePath
            << endl;
        return nullptr;
    }
    return vtk::Tools::readSurfaceFile(filePath);
}

} // End namespace Foam


// ************************************************************************* //
