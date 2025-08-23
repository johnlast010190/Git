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

Namespace
    Foam::functionObjects::runTimeVis::Utils

Description
    Collection of useful functions used by more than one class

SourceFiles
    Utils.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "primitives/ranges/MinMax/MinMax.H"
#include "vtkDataSet.h"
#include "vtkDataArray.h"
#include "vtkSmartPointer.h"
#include "vtkArrayCalculator.h"
#include <vector>

namespace Foam::functionObjects::runTimeVis::Utils {
    vtkSmartPointer<vtkArrayCalculator> getMagnitudeCalculator
    (
                const word& vector, const word& resultArrayName, bool isPointData
    );
    std::vector<Foam::scalar> linspace(Foam::scalar start, Foam::scalar end, label n);
    void getDataSetValue(vtkDataArray* dataArray, label id, label& size, scalar* array);
    void divideArray(const scalar* input, label size, scalar denominator, scalar* output);
    string getSubdictName(const dictionary& subdict);

    bool isNumericString(const std::string& s);

    void copyFieldData(vtkDataSet* source, vtkDataSet* destination);

    vtkSmartPointer<vtkPolyData> mergePolyData(vtkPolyData* polyData);
    vtkSmartPointer<vtkPolyData> getSurfacePolyDataFrom(vtkDataSet* dataSet);
    scalar calculateApproximateMaxElementSize(vtkDataSet* dataSet);

    bool hasAnyGhostPoint(vtkDataSet* dataSet);
    bool hasAnyGhostCell(vtkDataSet* dataSet);

    struct ArrayCorrespondence
    {
        vtkDataArray *baseArray;
        vtkDataArray *boundaryArray;

        ArrayCorrespondence(vtkDataArray *baseArray, vtkDataArray *boundaryArray)
            : baseArray(baseArray), boundaryArray(boundaryArray)
        {};
    };

    std::vector<ArrayCorrespondence> determineArrayCorrespondence(
        vtkDataSetAttributes *dataset,
        vtkDataSetAttributes *boundary,
        const std::vector<std::string>& excludedArrays
    );

    bool isRangeValid(scalarMinMax range);

    template<class T>
    bool contains(const std::vector<T>& vec, const T& key)
    {
        return std::find(vec.begin(), vec.end(), key) =! vec.end();
    }
} // End namespace Foam