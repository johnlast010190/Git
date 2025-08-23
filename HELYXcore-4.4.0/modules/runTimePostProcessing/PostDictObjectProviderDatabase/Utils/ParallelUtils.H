/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : Dev
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
    (c) 2023 Engys Ltd.

Namespace
    Foam::functionObjects::runTimeVis::ParallelUtils

Description
    Collection of useful parallel handling functions used by more than one class

SourceFiles
    ParallelUtils.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "primitives/ranges/MinMax/MinMax.H"
#include "vtkDataSet.h"
#include "vtkDataArray.h"
#include "vtkSmartPointer.h"
#include "vtkMultiProcessController.h"
#include "vtkArrayCalculator.h"
#include <vector>

namespace Foam::functionObjects::runTimeVis::ParallelUtils {
    void tradeValueWithProc(const scalar* input, scalar* output, label length, label procId);
    label getMaxLabelFromAllProcs(label value, vtkMultiProcessController* controller);
    void debugInfoParallel(const std::string& s);
    void debugInfoAllProcs(const string& text, scalar value);
    void debugInfoAllProcs(const string& text, const string& value);
    List<List<scalar>> gatherScalarList(const List<scalar>& values);
    List<List<scalar>> allGatherScalarList(const List<scalar>& values);
    List<List<scalar>> scatterScalarListList(const List<List<scalar>> &values);
    bool isRunningInParallel();
    bool isMaster();
    int localProcessId();
} // End namespace Foam