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
    meshSurfaceMapper2D

Description
    Maps vertices of a 2D onto the nearest points on the geometry

SourceFiles
    meshSurfaceMapper2D.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceMapper2D_H
#define meshSurfaceMapper2D_H

#include "primitives/ints/lists/labelList.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "utilities/containers/Lists/labelLongList.H"
#include "utilities/surfaceTools/meshSurfaceMapper/parMapperHelper.H"
#include "meshes/boundBox/boundBox.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class polyMeshGen;
class polyMeshGen2DEngine;
class meshSurfaceEngine;
class meshSurfacePartitioner;
class triSurfacePartitioner;

/*---------------------------------------------------------------------------*\
                    Class meshSurfaceMapper2D Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceMapper2D
{
    // Private data

        //- Reference to mesh surface
        const meshSurfaceEngine& surfaceEngine_;

        //- Reference to the octree
        const meshOctree& meshOctree_;

        //- Mesh 2D engine
        mutable polyMeshGen2DEngine* mesh2DEnginePtr_;

        //- TriSurface partitioner
        mutable triSurfacePartitioner* surfPartitionerPtr_;

        //- Mesh surface partitioner
        mutable meshSurfacePartitioner* meshPartitionerPtr_;

        //- Labels of edges which are at the boundary which shall be mapped
        //  onto the surface
        labelLongList activeBoundaryEdges_;


    // Private member functions

        //- Find active boundary edges
        void findActiveBoundaryEdges();

        //- Create a 2D engine
        void create2DEngine() const;

        //- Create and return mesh 2D engine
        inline const polyMeshGen2DEngine& mesh2DEngine() const
        {
            if (!mesh2DEnginePtr_)
                create2DEngine();

            return *mesh2DEnginePtr_;
        }

        //- Create surface partitioner
        void createTriSurfacePartitioner() const;

        //- Create and reurn surface partitioner
        inline const triSurfacePartitioner& surfacePartitioner() const
        {
            if (!surfPartitionerPtr_)
                createTriSurfacePartitioner();

            return *surfPartitionerPtr_;
        }

        //- Create mesh surface partitioner
        void createMeshSurfacePartitioner() const;

        //- Create and return mesh surface partitioner
        inline const meshSurfacePartitioner& meshPartitioner() const
        {
            if (!meshPartitionerPtr_)
                createMeshSurfacePartitioner();

            return *meshPartitionerPtr_;
        }

        //- Delete dynamically allocated data
        void clearOut();


    // Private member functions needed for parallel execution

        //- Find mapping distance for corner points
        void findMappingDistance
        (
            const labelLongList&,
            std::map<label, scalar>&
        ) const;

        //- Map to the location with the smallest distance
        void mapToSmallestDistance(LongList<parMapperHelper>&);


public:

    // Constructors

        //- Construct from meshSurfaceEngine and octree
        meshSurfaceMapper2D(const meshSurfaceEngine&, const meshOctree&);

        //- Disallow default bitwise copy construction
        meshSurfaceMapper2D(const meshSurfaceMapper2D&) = delete;


    //- Destructor
    ~meshSurfaceMapper2D();


    // Member Functions

        //- Adjust z coordinates of the mesh to the ones in the surface mesh
        void adjustZCoordinates();

        //- Projects surface vertices onto their nearest location
        //  on the surface mesh
        void mapVerticesOntoSurface();

        //- Projects selected edge to their nearest location
        //  on the surface mesh
        void mapVerticesOntoSurface(const labelLongList& edgesToMap);

        //- Projects edges with corners onto their nearest location
        //  on the surface mesh
        void mapCorners();

        //- Projects selected edges with corners onto their nearest
        //  locations on the surface mesh
        void mapCorners(const labelLongList& edgesToMap);

        //- Projects surface edges onto the surface with respect
        //  to the surface patch they belong to. Corners are respected
        void mapVerticesOntoSurfacePatches();

        //- Projects selected surface edges onto the surface with respect
        //  to the surface patch they belong to. Corners are respected
        void mapVerticesOntoSurfacePatches(const labelLongList& edgesToMap);

        //- A combination of mapping and smoothing intended for better
        //  feature capturing
        void preMapVertices(const label nIterations = 2);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceMapper2D&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
