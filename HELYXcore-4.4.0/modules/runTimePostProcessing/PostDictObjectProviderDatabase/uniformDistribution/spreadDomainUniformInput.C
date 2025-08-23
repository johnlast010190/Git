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
    (c) 2020-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "spreadDomainUniformInput.H"

#include "Utils/boundsUtils.H"
#include "Utils/ParallelUtils.H"

namespace Foam::functionObjects::runTimeVis
{

SpreadDomainUniformInput::SpreadDomainUniformInput
(
        vtkDataSet* inputData,
        label maximumNumberOfPoints,
        int seed
) : BoundsBasedUniformInput(inputData, maximumNumberOfPoints, seed)
{
    SpreadDomainUniformInput::computeBounds(bounds_);
    compliancyDistanceSquared_ = calculateCompliancyDistanceSquared(bounds_);
}

void SpreadDomainUniformInput::computeBounds(scalar bounds[6])
{
    boundsUtils::getFullPolyBounds(inputData_, bounds);
}

vtkSmartPointer<vtkIdTypeArray> SpreadDomainUniformInput::filterPointsInsideDomain
(
        const UnfilteredPointProperties& unfilteredPointProperties
) {
    vtkSmartPointer<vtkIdTypeArray> filteredPointIds = vtkSmartPointer<vtkIdTypeArray>::New();
    if (ParallelUtils::isRunningInParallel()) {

        // We need to filter out points that belong to other processors
        vtkMultiProcessController* controller = vtkMultiProcessController::GetGlobalController();
        int nProcs = controller->GetNumberOfProcesses();
        int thisProc = controller->GetLocalProcessId();

        List<List<scalar>> procPointDistances = ParallelUtils::allGatherScalarList(unfilteredPointProperties.distances);

        int minProc;
        for (label i = 0; i < unfilteredPointProperties.distances.size(); i++) {
            scalar minDistance = compliancyDistanceSquared_;
            minProc = -1;
            for (int j = 0; j < nProcs; j++) {
                // Compare the distance from the original point to the closest point in each processor domain
                if (procPointDistances[j][i] <= minDistance) {
                    minDistance = procPointDistances[j][i];
                    minProc = j;
                }
            }

            // If this processor point has the closest distance to the original point, and this distance is smaller than the compliance distance, then add it
            if (minProc == thisProc) {
                filteredPointIds->InsertNextValue(unfilteredPointProperties.ids->GetValue(i));
            }
        }
    } else {
        // If it's a serial case we just need to make sure that all points are within the compliance distance
        filteredPointIds = filterCompliantPoints(unfilteredPointProperties);
    }
    return filteredPointIds;
}

} // End namespace runTimeVis
// End namespace functionObjects
// End namespace Foam
