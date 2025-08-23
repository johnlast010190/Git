/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
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
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

Description
    Update the polyMesh corresponding to the given map.

\*---------------------------------------------------------------------------*/

#include "meshes/polyMesh/polyMesh.H"
#include "meshes/polyMesh/polyTopoChangeMap/polyTopoChangeMap.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/globalMeshData/globalMeshData.H"
#include "meshes/pointMesh/pointMesh.H"
#include "algorithms/indexedOctree/indexedOctree.H"
#include "algorithms/indexedOctree/treeDataCell.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyMesh::topoChange(const polyTopoChangeMap& map)
{
    if (debug)
    {
        InfoInFunction
            << "Updating addressing and (optional) pointMesh/pointFields"
            << endl;
    }

    // Update boundaryMesh (note that patches themselves already ok)
    boundary_.topoChange();

    // Update zones
    pointZones_.clearAddressing();
    faceZones_.clearAddressing();
    cellZones_.clearAddressing();

    // Remove the stored tet base points
    tetBasePtIsPtr_.clear();
    // Remove the cell tree
    cellTreePtr_.clear();

    // Update parallel data
    if (globalMeshDataPtr_.valid())
    {
        globalMeshDataPtr_->topoChange();
    }

    setInstance(time().timeName());

    // Map the old motion points if present
    if (oldPointsPtr_.valid())
    {
        // Make a copy of the original points
        pointField oldMotionPoints = oldPointsPtr_();

        pointField& newMotionPoints = oldPointsPtr_();

        // Resize the list to new size
        newMotionPoints.setSize(points_.size());

        // Map the list
        if (map.hasMotionPoints())
        {
            newMotionPoints.map(oldMotionPoints, map.pointMap());

            // Any points created out-of-nothing get set to the current
            // coordinate for lack of anything better.
            forAll(map.pointMap(), newPointi)
            {
                if (map.pointMap()[newPointi] == -1)
                {
                    newMotionPoints[newPointi] = points_[newPointi];
                }
            }
        }
        else
        {
            const labelList& pointMap = map.pointMap();
            const labelList& revPointMap = map.reversePointMap();

            forAll(pointMap, newPointi)
            {
                label oldPointi = pointMap[newPointi];
                if (oldPointi >= 0)
                {
                    if (revPointMap[oldPointi] == newPointi) // master point
                    {
                        newMotionPoints[newPointi] = oldMotionPoints[oldPointi];
                    }
                    else
                    {
                        newMotionPoints[newPointi] = points_[newPointi];
                    }
                }
                else
                {
                    newMotionPoints[newPointi] = points_[newPointi];
                }
            }
        }
    }

    meshObject::topoChange<polyMesh>(*this, map);
    meshObject::topoChange<pointMesh>(*this, map);

    // Reset valid directions (could change by faces put into empty patches)
    geometricD_ = Zero;
    solutionD_ = Zero;
}


void Foam::polyMesh::mapMesh(const polyMeshMap& map)
{
    meshObject::mapMesh<polyMesh>(*this, map);
    meshObject::mapMesh<pointMesh>(*this, map);
}


void Foam::polyMesh::distribute(const polyDistributionMap& map)
{
    meshObject::distribute<polyMesh>(*this, map);
    meshObject::distribute<pointMesh>(*this, map);
}


void Foam::polyMesh::updateGIB()
{
    calcDirections();
    boundary_.updateGIB();

    // Flags faceZones and cellZones files as being changed,
    // so they can be written down in the subsequent write time.
    faceZones_.instance() = time().timeName();
    cellZones_.instance() = time().timeName();

    faceZones_.writeOpt() = IOobject::AUTO_WRITE;
    cellZones_.writeOpt() = IOobject::AUTO_WRITE;
}


// ************************************************************************* //
