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
    meshOctreeCube

Description
    A cube stores information needed for mesh generation

SourceFiles
    meshOctreeCube.C

\*---------------------------------------------------------------------------*/

#ifndef meshOctreeCube_H
#define meshOctreeCube_H

#include "utilities/containers/DynList/DynList.H"
#include "utilities/containers/LongList/LongList.H"
#include "primitives/bools/lists/boolList.H"
#include "utilities/octrees/meshOctree/meshOctreeCube/meshOctreeCubeBasic.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Ostream;
class triSurf;
class VRWGraph;
class meshOctreeSlot;

/*---------------------------------------------------------------------------*\
                    Class meshOctreeCube Declaration
\*---------------------------------------------------------------------------*/

class meshOctreeCube
:
    public meshOctreeCubeBasic
{
    // Private data

        //- Pointer the slot containing this cube
        meshOctreeSlot* activeSlotPtr_;

        //- Pointer to the first child element
        meshOctreeCube** subCubesPtr_;

        //- Position of the cube in the list of leaves
        mutable label cubeLabel_;

        //- Label of the row which contains elements contained in the cube
        label containedElementsLabel_;

        //- Labels of contained surface edges
        label containedEdgesLabel_;


    // Static data

        //- This data is needed to map the the current refinement pattern
        //  to the Hilbert space-filling curve in order to reduce the bandwidth
        //  in the list of leaves. The procedure is performed using the table
        //  published in (Campbell etal: Dynamic Octree Load Balancing Using
        //  Space-Filling Curves) which is adjusted to the current ref pattern
        static const label hOrder_[24][8];

        //- Orientation of the Hilbert curve in the refined cubes of a cube
        //  with the ordering in the above table
        static const label hOrient_[24][8];


    // Private member functions

        //- Find edges contained in the cube
        void findContainedEdges
        (
            const triSurf&,
            const boundBox&
        );


public:

    // Constructors

        //- Default constructor
        inline meshOctreeCube();

        //- Copy construct
        inline meshOctreeCube(const meshOctreeCube&);

        //- Construct from coordinates
        meshOctreeCube(const meshOctreeCubeCoordinates&);

        //- Construct from coordinates and the number of surface triangles
        //- this constructor is used for the initial cube, only
        meshOctreeCube
        (
            const meshOctreeCubeCoordinates& cc,
            const label nElmts,
            meshOctreeSlot* slotPtr
        );

        //- Disallow default bitwise copy construction
        meshOctreeCube(const meshOctreeCubeBasic&) = delete;


    //- Destructor
    ~meshOctreeCube();


    // Member Functions

        //- Refine cube in two directions, it is used for generating quadtrees
        void refineCube2D
        (
            const triSurf&,
            const boundBox&,
            meshOctreeSlot* slotPtr = nullptr
        );

        //- Subdivide the octree cube
        void refineCube
        (
            const triSurf&,
            const boundBox&,
            meshOctreeSlot* slotPtr = nullptr
        );

        //- Create missing child cubes and refine them until the required
        //  cube is reached (this function is used for parallel octree creation)
        void refineMissingCube
        (
            const triSurf&,
            const boundBox&,
            const label scI,
            meshOctreeSlot* slotPtr = nullptr
        );

        //- Create missing child cubes and refine them until the required
        //  cube is reached (this function is used for parallel octree creation)
        //  contained elements and edges are provided as argument
        void refineMissingCube
        (
            const label scI,
            const label elementsRowI = -1,
            const label edgesRowI = -1,
            meshOctreeSlot* slotPtr = nullptr
        );

        //- Return the pointer to the slot containing the cube
        inline const meshOctreeSlot* slotPtr() const;

        //- Position of the cube in the list of leaves
        inline label cubeLabel() const;

        //- Check if the cube is a leaf
        inline bool isLeaf() const;

        //- Return a pointer to a child cubes at given position
        inline meshOctreeCube* subCube(const label) const;

        //- Return the pointers to the sons
        FixedList<meshOctreeCube*, 8> subCubes() const;

        //- Check if this box has some contained triangles
        bool hasContainedTriangles
        (
            const triSurf&,
            const boundBox&,
            const VRWGraph& containedElements
        ) const;

        //- Return true if the box contains some triangles
        inline bool hasContainedElements() const;

        //- Returns the list of contained surface elements
        //  this list is only allocated for cubes containing elements
        inline label containedElements() const;

        //- Returns true if the box contains some edges
        inline bool hasContainedEdges() const;

        //- Returns the list of contained surface edges
        //  this list is only allocated for cubes containing edges
        inline label containedEdges() const;


    // Recursive functions

        //- Leaves contained in the given box
        void leavesInBox
        (
            const boundBox& rootBox,
            const boundBox& searchingBox,
            DynList<const meshOctreeCube*, 256>&
        ) const;

        //- Find leaves within a sphere
        void leavesInSphere
        (
            const boundBox& rootBox,
            const point&,
            const scalar,
            DynList<label>&
        ) const;

        //- Mark leaves within a sphere
        void markLeavesInSphere
        (
            const boundBox& rootBox,
            const point&,
            const scalar,
            labelList&,
            bool&
        ) const;

        //- Find leaves for a given cube
        void findLeaves
        (
            LongList<meshOctreeCube*>& leaves
        ) const;

        //- Find coordinates of cubes which are located on other processors
        void findCoordinatesOfMissingCubes
        (
            LongList<meshOctreeCubeCoordinates>& coordinates
        ) const;

        //- Count number of originating from this cube
        void countChildCubes(label& nCubes) const;

        //- Delete boxes which are not local to the given processor
        bool purgeProcessorCubes(const short procNo);


    // Member operators

        //- Assignment
        inline void operator=(const meshOctreeCube&);
        friend Ostream& operator<<(Ostream&, const meshOctreeCube&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "utilities/octrees/meshOctree/meshOctreeCube/meshOctreeCubeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
