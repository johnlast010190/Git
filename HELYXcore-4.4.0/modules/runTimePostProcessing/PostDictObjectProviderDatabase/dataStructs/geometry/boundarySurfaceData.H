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
    Foam::functionObjects::runTimeVis::BoundarySurfaceData

Description
    Contains the dict data necessary for the boundaries of a user defined block mesh

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "hash/hasher.H"
#include "meshes/primitiveShapes/point/point.H"

#include "dataStructs/geometry/blockMesh/blockMeshBlock.H"
#include "dataStructs/geometry/blockMesh/blockMeshPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

struct BoundarySurfaceData {
    static const constexpr int NUM_POINTS = 4;
    point points[NUM_POINTS];
    label elements[2] = {0, 0};

    BoundarySurfaceData(const BlockMeshPatch& patch, const BlockMeshBlock& block, const List<point>& vertices)
    {
        if (patch.type == BlockMeshPatchType::WALL && block.type == BlockMeshBlockType::HEX)
        {
            if (patch.faces.size() == 1 && patch.faces[0].size() == NUM_POINTS)
            {
                for (int i = 0; i < NUM_POINTS; i++)
                {
                    this->points[i] = vertices[patch.faces[0][i]];
                }
                point axis1 = points[1] - points[0];
                axis1.normalise();
                point axis2 = points[3] - points[0];
                axis2.normalise();
                point elementsPoint(block.elements[0], block.elements[1], block.elements[2]);
                elements[0] = label(round(dot(elementsPoint, axis1)));
                elements[1] = label(round(dot(elementsPoint, axis2)));
            }
        }
        else
        {
            WarningInFunction << "Unknown block mesh boundary surface type combination for " << patch.name << endl;
        }
    }

    bool operator==(const BoundarySurfaceData& other) const
    {
        for (int i = 0; i < NUM_POINTS; i++)
        {
            if (points[i] != other.points[i])
            {
                return false;
            }
        }
        for (label i = 0; i < 2; i ++)
        {
            if (elements[i] != other.elements[i])
            {
                return false;
            }
        }
        return true;
    }

    [[nodiscard]] size_t computeAndReturnHash() const
    {
        size_t hash = 0;
        for (const point& point : points)
        {
            hasher::hash_combine(hash, point);
        }
        for (int element : elements)
        {
            hasher::hash_combine(hash, element);
        }
        return hash;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
