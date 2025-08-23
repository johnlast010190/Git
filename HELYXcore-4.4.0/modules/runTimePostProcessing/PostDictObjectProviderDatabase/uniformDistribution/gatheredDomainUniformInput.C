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

#include "gatheredDomainUniformInput.H"

namespace Foam::functionObjects::runTimeVis
{

GatheredDomainUniformInput::GatheredDomainUniformInput
(
        vtkDataSet* inputData,
        label maximumNumberOfPoints,
        int seed
) : BoundsBasedUniformInput(inputData, maximumNumberOfPoints, seed)
{
    GatheredDomainUniformInput::computeBounds(bounds_);
    compliancyDistanceSquared_ = calculateCompliancyDistanceSquared(bounds_);
}

void GatheredDomainUniformInput::computeBounds(scalar bounds[6])
{
#if defined(HELYX_SP)
    double dBounds[6];
    inputData_->GetBounds(dBounds);
    for (int i = 0; i < 6; i++) bounds[i] = static_cast<scalar>(dBounds[i]);
#else
    inputData_->GetBounds(bounds);
#endif
}

vtkSmartPointer<vtkIdTypeArray> GatheredDomainUniformInput::filterPointsInsideDomain
(
        const UnfilteredPointProperties& unfilteredPointProperties
)
{
    return filterCompliantPoints(unfilteredPointProperties);
}

} // End namespace runTimeVis
// End namespace functionObjects
// End namespace Foam
