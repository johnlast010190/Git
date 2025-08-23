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
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::ThresholdDataSetProvider

Description
    Provider for the threshold 3d objects, which computes their vtk dataset

SourceFiles
    thresholdDataSetProvider.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "itemDataSetProviders/itemDataSetProvider.H"

#include "hash/hasher.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class engysThreshold;

namespace Foam::functionObjects::runTimeVis
{

struct ThresholdObjectData
{
    foamField field;
    bool allPointsCriterion;
    scalar minThreshold;
    scalar maxThreshold;
    bool invert;

    bool operator==(const ThresholdObjectData& other) const
    {
        return           field == other.field &&
            allPointsCriterion == other.allPointsCriterion &&
                  minThreshold == other.minThreshold &&
                  maxThreshold == other.maxThreshold &&
                        invert == other.invert;
    }

    [[nodiscard]] size_t computeAndReturnHash() const
    {
        size_t hash = 0;
        hasher::hash_combine(hash, field);
        hasher::hash_combine(hash, allPointsCriterion);
        hasher::hash_combine(hash, minThreshold);
        hasher::hash_combine(hash, maxThreshold);
        hasher::hash_combine(hash, invert);
        return hash;
    }
};

/*---------------------------------------------------------------------------*\
                 Class TransformDataSetProvider Declaration
\*---------------------------------------------------------------------------*/

class ThresholdDataSetProvider : public ItemDataSetProvider
{
private:

protected:

        vtkSmartPointer<engysThreshold> thresholdFilter_;

        void update(scalar currentTime) override;

public:
        ThresholdDataSetProvider(const std::string& name, const ThresholdObjectData& dictData);
        ~ThresholdDataSetProvider() override = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
