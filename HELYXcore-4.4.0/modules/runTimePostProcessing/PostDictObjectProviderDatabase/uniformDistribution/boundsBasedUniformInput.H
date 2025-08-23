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

Class
    Foam::functionObjects::runTimeVis::BoundsBasedUniformInput

Description
    Base class for calculating the bounds based uniform input seed type for glyphs and streamlines.
    Use one of SpreadDomainUniformInput or GatheredDomainUniformInput instead

SourceFiles
    boundsBasedUniformInput.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimeVis_boundsBasedUniformInput_H
#define functionObjects_runTimeVis_boundsBasedUniformInput_H

#include "meshes/boundBox/boundBox.H"
#include "vtkPolyData.h"

#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class vtkDataObject;
class vtkDataSet;

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
                   Class BoundsBasedUniformInput Declaration
\*---------------------------------------------------------------------------*/
class BoundsBasedUniformInput
{
protected:
    vtkSmartPointer<vtkDataSet> inputData_;
    label maximumNumberOfPoints_;
    int seed_;
    scalar bounds_[6];

    scalar compliancyDistanceSquared_;

    struct UnfilteredPointProperties {
        vtkSmartPointer<vtkIdTypeArray> ids;
        DynamicList<scalar> distances;
    };

public:
    BoundsBasedUniformInput(
            vtkDataSet* inputData,
            label maximumNumberOfPoints,
            int seed
    );

    // Returns a number of points of the input data uniformly distributed over its bounds, using a given seed
    vtkSmartPointer<vtkPolyData> compute();

protected:
    virtual void computeBounds(scalar bounds[6]) = 0;

    [[nodiscard]] virtual vtkSmartPointer<vtkIdTypeArray> filterPointsInsideDomain
    (
        const UnfilteredPointProperties& unfilteredPointProperties
    ) = 0;

    [[nodiscard]] vtkSmartPointer<vtkIdTypeArray> filterCompliantPoints
    (
        const UnfilteredPointProperties& unfilteredPointProperties
    ) const;

    [[nodiscard]] scalar calculateCompliancyDistanceSquared(const scalar bounds[6]) const;

    virtual ~BoundsBasedUniformInput() = 0;

private:
    List<Vector<scalar>> generateUniformPointDistribution();

    vtkSmartPointer<vtkIdTypeArray> findPointsInsideDomain
    (
        const List<Vector<scalar>>& customPointDistribution
    );

    vtkSmartPointer<vtkPolyData> extractInputPointsFromPointIds
    (
            const vtkSmartPointer<vtkIdTypeArray> &filteredPointIds
    );

    UnfilteredPointProperties findClosestPointsAndDistancesToThem
    (
        const List<Vector<scalar>>& customPointDistribution
    );

};

} // End namespace

#endif
