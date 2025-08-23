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
    (c) Creative Fields, Ltd.
    (c) 2024 Engys Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Class
    meshSurfaceMapper

Description
    Maps vertices to the nearest point on the geometry surface

SourceFiles
    meshSurfaceMapper.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceMapper_H
#define meshSurfaceMapper_H

#include "primitives/ints/lists/labelList.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "utilities/containers/Lists/labelLongList.H"
#include "utilities/surfaceTools/meshSurfaceMapper/parMapperHelper.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class meshSurfaceEngine;
class meshSurfacePartitioner;
class triSurfacePartitioner;

/*---------------------------------------------------------------------------*\
                    Class meshSurfaceMapper Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceMapper
{
    // Private data

        //- Mesh surface
        const meshSurfaceEngine& surfaceEngine_;

        //- Reference to the octree
        const meshOctree& meshOctree_;

        //- Mesh surface partitioner
        mutable const meshSurfacePartitioner* surfaceEnginePartitionerPtr_;
        const bool deletePartitioner_;

        //- TriSurface partitioner
        mutable triSurfacePartitioner* surfPartitionerPtr_;


    // Private member functions

        //- Create and return mesh surface partitioner
        void createMeshSurfacePartitioner() const;
        inline const meshSurfacePartitioner& meshPartitioner() const
        {
            if (!surfaceEnginePartitionerPtr_)
                createMeshSurfacePartitioner();

            return *surfaceEnginePartitionerPtr_;
        }

        //- Create and return triSurfacePartitioner
        void createTriSurfacePartitioner() const;
        inline const triSurfacePartitioner& surfacePartitioner() const
        {
            if (!surfPartitionerPtr_)
                createTriSurfacePartitioner();

            return *surfPartitionerPtr_;
        }

        //- Delete surfaceEnginePartitionerPtr_ and surfPartitionerPtr_
        void clearOut();


    // Private member functions

        //- Map corner nodes to the boundary
        void mapCorners(const labelLongList& nodesToMap);

        //- Find mapping distance for selected points
        void findMappingDistance
        (
            const labelLongList& nodesToMap,
            scalarList & mappingDistance
        ) const;

        //- Calculate face metric
        scalar faceMetricInPatch(const label bfI, const label patchI) const;


    // Private member functions needed for parallel execution

        //- Check if nodes at parallel boundaries are selected at all processors
        void selectNodesAtParallelBnd(const labelLongList&);

        //- Map to the smallest distance
        //  makes sense for parallel calculations, only
        void mapToSmallestDistance(LongList<parMapperHelper>&);


public:

    // Constructors

        //- Construct from meshSurfaceEngine and octree
        meshSurfaceMapper(const meshSurfaceEngine&, const meshOctree&);

        //- Construct from meshSurfacePartitioner and octree
        meshSurfaceMapper(const meshSurfacePartitioner&, const meshOctree&);

        //- Disallow default bitwise copy construction
        meshSurfaceMapper(const meshSurfaceMapper&) = delete;


    //- Destructor
    ~meshSurfaceMapper();


    // Member Functions

        //- Map the given surface node on the selected patch
        //  this does not make sense for vertices at parallel boundaries
        void mapNodeToPatch(const label bpI, const label patchI = -1);

        //- Projects surface vertices onto their nearest location
        //  on the surface mesh
        void mapVerticesOntoSurface();

        //- Projects selected surface vertices to their nearest location
        //  on the surface mesh
        void mapVerticesOntoSurface(const labelLongList& nodesToMap);

        //- Projects corner and edge vertices onto their nearest location
        //  on the surface mesh
        void mapCornersAndEdges();

        //- Projects selected edge vertices onto their nearest
        //  locations on the surface mesh
        void mapEdgeNodes(const labelLongList& nodesToMap);

        //- Projects surface vertices onto the surface with respect
        //  to the surface patch they belong to. Edges and corner are respected
        void mapVerticesOntoSurfacePatches();

        //- Projects selected surface vertices onto the surface with respect
        //  to the surface patch they belong to. Edges and corner are respected
        void mapVerticesOntoSurfacePatches(const labelLongList& nodesToMap);

        //- A combination of mapping and smoothing intended for better
        //  feature capturing
        void preMapVertices(const label nIterations = 3);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceMapper&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
