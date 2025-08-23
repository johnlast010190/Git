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

#include "boundsBasedUniformInput.H"

#include "containers/Lists/List/List.H"

#include "vtkMinimalStandardRandomSequence.h"
#include "vtkOctreePointLocator.h"
#include "vtkIdTypeArray.h"
#include "vtkSelection.h"
#include "vtkSelectionNode.h"
#include "vtkExtractSelection.h"
#include "vtkGeometryFilter.h"

namespace Foam::functionObjects::runTimeVis
{

BoundsBasedUniformInput::BoundsBasedUniformInput
(
        vtkDataSet* inputData,
        label maximumNumberOfPoints,
        int seed
) :
    inputData_(inputData),
    maximumNumberOfPoints_(maximumNumberOfPoints),
    seed_(seed),
    bounds_(),
    compliancyDistanceSquared_(0)
{
}

BoundsBasedUniformInput::~BoundsBasedUniformInput() = default;

vtkSmartPointer<vtkPolyData> BoundsBasedUniformInput::compute()
{
    // Let every proc generate all points, based on a universal seed.  This
    // is less daft than it sounds - either every proc just waits there, and
    // we have to broadcast the results, or every proc just does the
    // computation
    List<Vector<scalar>> points = generateUniformPointDistribution();
    vtkSmartPointer<vtkIdTypeArray> filteredIds = findPointsInsideDomain(points);
    return extractInputPointsFromPointIds(filteredIds);
}

List<Vector<scalar>> BoundsBasedUniformInput::generateUniformPointDistribution()
{
    List<Vector<scalar>> points(maximumNumberOfPoints_);

    vtkNew<vtkMinimalStandardRandomSequence> randomGenerator;
    randomGenerator->SetSeed(seed_);

    for (auto & point : points)
    {
        randomGenerator->Next();
        point[0] = static_cast<scalar>(randomGenerator->GetRangeValue(bounds_[0], bounds_[1]));
        randomGenerator->Next();
        point[1] = static_cast<scalar>(randomGenerator->GetRangeValue(bounds_[2], bounds_[3]));
        randomGenerator->Next();
        point[2] = static_cast<scalar>(randomGenerator->GetRangeValue(bounds_[4], bounds_[5]));
    }

    return points;
}

scalar BoundsBasedUniformInput::calculateCompliancyDistanceSquared(const scalar bounds[6]) const
{
    Vector<scalar> lengthsWithSign;
    lengthsWithSign[0] = bounds[1] - bounds[0];
    lengthsWithSign[1] = bounds[3] - bounds[2];
    lengthsWithSign[2] = bounds[5] - bounds[4];

    label dimensions =
        (
            lengthsWithSign[0] != 0.0 &&
            lengthsWithSign[1] != 0.0 &&
            lengthsWithSign[2] != 0.0
        ) ? 3 : 2;
    // Feel like there's probably a build-in for this...
    scalar diagonalLength = mag(lengthsWithSign);
    auto volume = static_cast<scalar>(std::pow(diagonalLength, dimensions));

    scalar nearestPointRadius = 0.0001;
    if (volume > 0.0)
    {
        // numberOfPoints cannot be zero - rule enforced in GUI
        scalar volumePerGlyph = volume / static_cast<scalar>(maximumNumberOfPoints_);
        auto delta = static_cast<scalar>(pow(volumePerGlyph, 1.0 / dimensions));
        nearestPointRadius = delta / static_cast<scalar>(2.0);
    }
    return nearestPointRadius * nearestPointRadius;
}

vtkSmartPointer<vtkIdTypeArray> BoundsBasedUniformInput::findPointsInsideDomain
(
    const List<Vector<scalar>>& customPointDistribution
)
{
    UnfilteredPointProperties unfilteredPointProperties = findClosestPointsAndDistancesToThem(customPointDistribution);
    return filterPointsInsideDomain(unfilteredPointProperties);
}

BoundsBasedUniformInput::UnfilteredPointProperties BoundsBasedUniformInput::findClosestPointsAndDistancesToThem
(
    const List<Vector<scalar>>& customPointDistribution
)
{
    // At the moment, this does nothing!  Need to break down into different
    // bounding boxes
    UnfilteredPointProperties result;
    result.ids = vtkSmartPointer<vtkIdTypeArray>::New();

    if (inputData_->GetNumberOfPoints() > 0)
    {
        vtkSmartPointer<vtkOctreePointLocator> pointLocator =
        vtkSmartPointer<vtkOctreePointLocator>::New();
        pointLocator->Initialize();
        pointLocator->SetDataSet(inputData_);
        pointLocator->BuildLocator();

        for (auto foamPoint : customPointDistribution)
        {
            double point[3];
            for (auto j = 0; j < 3; j++)
            {
                point[j] = foamPoint[j];
            }

            // label closestPointId = pointLocator->FindClosestPoint(point);
            vtkIdType closestPointId = pointLocator->FindClosestPoint(point);
            result.ids->InsertNextValue(closestPointId);

            double closestPoint[3];
            inputData_->GetPoint(closestPointId, closestPoint);

            // Squared distance between points
            Vector<scalar> testPoint;
            for (auto j = 0; j < 3; j++)
            {
                testPoint[j] = static_cast<scalar>(closestPoint[j]);
            }

            result.distances.append(magSqr(foamPoint - testPoint));
        }
    } else {
        // Must add these distances because all processes need to have some distance values to avoid deadlocks
        for (label i = 0; i < customPointDistribution.size(); i++)
        {
            result.distances.append(compliancyDistanceSquared_ * 2 + 1);
        }
    }
    return result;
}

vtkSmartPointer<vtkIdTypeArray> BoundsBasedUniformInput::filterCompliantPoints
(
        const UnfilteredPointProperties& unfilteredPointProperties
) const
{
    vtkSmartPointer<vtkIdTypeArray> filteredPointIds =
        vtkSmartPointer<vtkIdTypeArray>::New();
    for (label i = 0; i < unfilteredPointProperties.distances.size(); i++)
    {
        if (unfilteredPointProperties.distances[i] <= compliancyDistanceSquared_)
        {
            filteredPointIds->InsertNextValue(unfilteredPointProperties.ids->GetValue(i));
        }
    }
    return filteredPointIds;
}

vtkSmartPointer<vtkPolyData> BoundsBasedUniformInput::extractInputPointsFromPointIds
(
        const vtkSmartPointer<vtkIdTypeArray> &filteredPointIds
)
{
    if (filteredPointIds->GetNumberOfComponents() > 0)
    {

        vtkNew<vtkSelectionNode> node;
        node->SetContentType(vtkSelectionNode::SelectionContent::INDICES);
        node->SetFieldType(vtkSelectionNode::SelectionField::POINT);
        node->SetSelectionList(filteredPointIds);

        vtkNew<vtkSelection> selection;
        selection->AddNode(node);

        vtkNew<vtkExtractSelection> filter;
        filter->SetInputData(0, inputData_);
        filter->SetInputData(1, selection);
        filter->Update();

        vtkNew<vtkGeometryFilter> geometryFilter;
        geometryFilter->SetInputData(filter->GetOutput());
        geometryFilter->Update();
        return geometryFilter->GetOutput();
    }
    else
    {
        return vtkSmartPointer<vtkPolyData>::New();
    }
}

} // End namespace runTimeVis
// End namespace functionObjects
// End namespace Foam
