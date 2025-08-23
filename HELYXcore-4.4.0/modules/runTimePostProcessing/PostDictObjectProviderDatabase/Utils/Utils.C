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
    (c) 2022-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "Utils.H"

#include "db/dictionary/dictionary.H"

#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkImplicitFunction.h"
#include "vtkCleanPolyData.h"
#include "vtkIncrementalOctreePointLocator.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDoubleArray.h"

vtkSmartPointer<vtkArrayCalculator> Foam::functionObjects::runTimeVis::Utils::getMagnitudeCalculator
(
    const word& vector, const word& resultArrayName, bool isPointData
)
{
    vtkSmartPointer<vtkArrayCalculator> calculator =
        vtkSmartPointer<vtkArrayCalculator>::New();
    if (isPointData)
    {
        calculator->SetAttributeTypeToPointData();
    }
    else
    {
        calculator->SetAttributeTypeToCellData();
    }
    calculator->AddVectorArrayName(vector.c_str(), 0, 1, 2);
    word temp = "mag(" + vector + ")";
    calculator->SetFunction(temp.c_str());
    calculator->SetResultArrayName(resultArrayName.c_str());
    return calculator;
}

std::vector<Foam::scalar> Foam::functionObjects::runTimeVis::Utils::linspace(Foam::scalar start, Foam::scalar end, label n) {
    if (n >= 2) {
        std::vector<Foam::scalar> result(n);
        Foam::scalar increment = (end - start) / static_cast<scalar>(n - 1);

        for (label i = 0; i < n; i++) {
            result.at(i) = start + increment * static_cast<scalar>(i);
        }
        return result;
    } else {
        std::vector<Foam::scalar> result(1);
        result.at(0) = (end + start) / 2;
        return result;
    }
}

void Foam::functionObjects::runTimeVis::Utils::getDataSetValue(vtkDataArray* dataArray, label id, label& size, scalar* array) {
    size = dataArray->GetNumberOfComponents();
#if defined(HELYX_SP)
    std::vector<double> dArray(size);
    dataArray->GetTuple(id, dArray.data());
    for (label i = 0; i < size; i++) array[i] = static_cast<scalar>(dArray[i]);
#else
    dataArray->GetTuple(id, array);
#endif
    switch (size) {
        case 1:
        case 3:
        case 6:
        case 9:
            return;
        default:
            size = 0;
            return;
    }
}

void Foam::functionObjects::runTimeVis::Utils::divideArray(const scalar* input, label size, scalar denominator, scalar* output) {
    for (label i = 0; i < size; i++) {
        if (denominator != 0) {
            output[i] = input[i] / denominator;
        } else {
            output[i] = 0;
        }
    }
}

Foam::string Foam::functionObjects::runTimeVis::Utils::getSubdictName(const dictionary& subdict)
{
    size_t parentNameLength = subdict.parent().empty() || subdict.parent().name().empty() ? 0 : subdict.parent().name().length();
    if (parentNameLength <= 0)
    {
        return subdict.name().name();
    }
    else
    {
        if (parentNameLength + 1 > subdict.name().size())
        {
            FatalError << "Trying to get subdict name from " << subdict.name()
                << " with parent" << subdict.parent().name() << endl
                << abort(FatalError);
        }
        return subdict.name().substr(parentNameLength+1);
    }
}

bool Foam::functionObjects::runTimeVis::Utils::isNumericString(const std::string& s)
{
    return !s.empty() && std::find_if(
        s.begin(),
        s.end(), [](unsigned char c) { return !std::isdigit(c); }
    ) == s.end();
}

void Foam::functionObjects::runTimeVis::Utils::copyFieldData(vtkDataSet* source, vtkDataSet* destination)
{
    vtkFieldData* sourceFieldData = source->GetFieldData();
    if (!sourceFieldData)
    {
        return;
    }

    if (!destination->GetFieldData())
    {
        vtkNew<vtkFieldData> newFieldData;
        destination->SetFieldData(newFieldData);
    }

    for (int i = 0; i < sourceFieldData->GetNumberOfArrays(); i++)
    {
        destination->GetFieldData()->AddArray(sourceFieldData->GetArray(i));
    }
}

std::vector<Foam::functionObjects::runTimeVis::Utils::ArrayCorrespondence>
Foam::functionObjects::runTimeVis::Utils::determineArrayCorrespondence(
    vtkDataSetAttributes *dataset,
    vtkDataSetAttributes *boundary,
    const std::vector<std::string>& excludedArrays
)
{
    vtkDataArray* boundaryGhostArray = boundary->GetArray(vtkCellData::GhostArrayName());
    std::vector<ArrayCorrespondence> correspondence;
    for (int boundaryArrayId = 0; boundaryArrayId < boundary->GetNumberOfArrays(); boundaryArrayId++)
    {
        vtkDataArray* boundaryArray = boundary->GetArray(boundaryArrayId);
        if (boundaryGhostArray == boundaryArray) continue;
        vtkDataArray* baseArray = dataset->GetArray(boundaryArray->GetName());
        if (baseArray &&
            std::find(excludedArrays.begin(), excludedArrays.end(), boundaryArray->GetName()) == excludedArrays.end())
        {
            correspondence.emplace_back(baseArray, boundaryArray);
        }
    }
    return correspondence;
}

bool Foam::functionObjects::runTimeVis::Utils::isRangeValid(scalarMinMax range)
{
    return range.valid() && !std::isnan(range.min()) && !std::isnan(range.max());
}

vtkSmartPointer<vtkPolyData> Foam::functionObjects::runTimeVis::Utils::mergePolyData(vtkPolyData* polyData)
{
    vtkNew<vtkIncrementalOctreePointLocator> locator;
    locator->SetTolerance(0);

    vtkNew<vtkCleanPolyData> mergeFilter;
    mergeFilter->SetInputData(polyData);
    mergeFilter->SetTolerance(0);
    mergeFilter->SetLocator(locator);
    mergeFilter->Update();
    return mergeFilter->GetOutput();
}

vtkSmartPointer<vtkPolyData> Foam::functionObjects::runTimeVis::Utils::getSurfacePolyDataFrom(vtkDataSet* dataSet)
{
    if (dataSet->IsA("vtkPolyData")) {
        vtkNew<vtkPolyData> shallowCopy;
        shallowCopy->ShallowCopy(dataSet);
        return shallowCopy;
    } else {
        vtkNew<vtkDataSetSurfaceFilter> surfaceFilter;
        surfaceFilter->SetInputData(dataSet);
        surfaceFilter->Update();
        return surfaceFilter->GetOutput();
    }
}
