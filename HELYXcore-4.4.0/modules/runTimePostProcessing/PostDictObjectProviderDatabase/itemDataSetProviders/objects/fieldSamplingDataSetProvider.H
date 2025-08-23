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
    (c) 2022 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::FieldSamplingDataSetProvider

Description
    Provider for the field sampling 3d objects, which computes their vtk dataset

SourceFiles
    fieldSamplingDataSetProvider.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "itemDataSetProviders/itemDataSetProvider.H"

#include "hash/hasher.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class engysResampleToImage;

namespace Foam::functionObjects::runTimeVis
{
struct FieldSamplingObjectData;

/*---------------------------------------------------------------------------*\
                 Class FieldSamplingDataSetProvider Declaration
\*---------------------------------------------------------------------------*/

class FieldSamplingDataSetProvider : public ItemDataSetProvider
{
private:

protected:

    vtkSmartPointer<engysResampleToImage> fieldSamplingFilter_;

    void update(scalar currentTime) override;

    void addExtraRequirementsDueToUpstreamRequirements(
        ItemRequirements &ownRequirements,
        const ItemRequirements &upstreamRequirements
    ) const override;

public:
    FieldSamplingDataSetProvider(const std::string& name, const FieldSamplingObjectData& dictData);
    ~FieldSamplingDataSetProvider() override = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
