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
    Foam::functionObjects::runTimeVis::GatheredDomainUniformInput

Description
    Class for calculating the bounds based uniform input seed type for streamlines,
     where each process has the whole domain.

SourceFiles
    gatheredDomainUniformInput.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "uniformDistribution/boundsBasedUniformInput.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
                 Class GatheredDomainUniformInput Declaration
\*---------------------------------------------------------------------------*/
class GatheredDomainUniformInput : public BoundsBasedUniformInput
{
public:
    GatheredDomainUniformInput(
            vtkDataSet* inputData,
            label maximumNumberOfPoints,
            int seed
    );

private:
    void computeBounds(scalar bounds[6]) override;

    [[nodiscard]] vtkSmartPointer<vtkIdTypeArray> filterPointsInsideDomain
    (
        const UnfilteredPointProperties& unfilteredPointProperties
    ) override;

};

} // End namespace runTimeVis
// End namespace functionObjects
// End namespace Foam
