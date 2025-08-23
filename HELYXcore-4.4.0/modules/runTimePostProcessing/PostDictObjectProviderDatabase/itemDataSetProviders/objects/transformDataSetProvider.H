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
    Foam::functionObjects::runTimeVis::TransformDataSetProvider

Description
    Provider for the transform 3d objects, which computes their vtk dataset

SourceFiles
    transformDataSetProvider.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "itemDataSetProviders/itemDataSetProvider.H"
#include "meshes/primitiveShapes/point/point.H"
#include "vector"
#include "types/transformationType.H"
#include "hash/hasher.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class vtkTransformFilter;
class vtkAppendFilter;

namespace Foam::functionObjects::runTimeVis
{

struct SingleTransformationObjectData
{
    TransformationType type;

    bool includeSource = false;
    point transformPoint;
    vector transformVector;
    scalar transformOffset = 0;
    label transformCount = 0;

    bool operator==(const SingleTransformationObjectData& other) const
    {
        return type == other.type &&
               includeSource == other.includeSource &&
               transformPoint == other.transformPoint &&
               transformVector == other.transformVector &&
               transformOffset == other.transformOffset &&
               transformCount == other.transformCount;
    }

    void computeAndAddToHash(size_t& hash) const
    {
        hasher::hash_combine(hash, type.getValue());
        hasher::hash_combine(hash, transformCount);
        hasher::hash_combine(hash, transformOffset);
        hasher::hash_combine(hash, transformPoint);
        hasher::hash_combine(hash, transformVector);
    }
};

struct TransformObjectData
{
    List<SingleTransformationObjectData> transformations;

    bool operator==(const TransformObjectData& other) const
    {
        if (transformations.size() != other.transformations.size())
        {
            return false;
        }
        // Here the order of the operations matters
        for (label i = 0; i < transformations.size(); i++)
        {
            if (!(transformations[i] == other.transformations[i]))
            {
                return false;
            }
        }
        return true;
    }

    [[nodiscard]] size_t computeAndReturnHash() const
    {
        size_t hash = 0;
        // Here the order of the transformations matters
        for (const SingleTransformationObjectData& transformation : transformations)
        {
            transformation.computeAndAddToHash(hash);
        }
        return hash;
    }
};

/*---------------------------------------------------------------------------*\
                   Class ItemDataSetProvider Declaration
\*---------------------------------------------------------------------------*/

class TransformDataSetProvider : public ItemDataSetProvider
{
private:

protected:

        std::vector<vtkSmartPointer<vtkTransformFilter>> transformFilters_;
        vtkSmartPointer<vtkAppendFilter> appendFilter_;

        void update(scalar currentTime) override;

public:
        TransformDataSetProvider(const std::string& name, const TransformObjectData& dictData);
        ~TransformDataSetProvider() override = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
