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
    meshOctreeAddressing

Description
    Automatic octree for mesh generation. This class checks if the current
    octree resolution is enough to resolve all surface features

SourceFiles
    meshOctreeAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef meshOctreeAddressing_H
#define meshOctreeAddressing_H

#include "utilities/containers/FRWGraph/FRWGraph.H"
#include "utilities/containers/VRWGraph/VRWGraph.H"
#include "utilities/octrees/meshOctree/meshOctreeCube/meshOctreeCubeCoordinates.H"
#include "meshes/meshShapes/face/face.H"
#include "meshes/meshShapes/edge/edge.H"
#include "utilities/containers/DynList/DynList.H"
#include "containers/HashTables/Map/Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class IOdictionary;

/*---------------------------------------------------------------------------*\
                    Class meshOctreeAddressing Declaration
\*---------------------------------------------------------------------------*/

class meshOctreeAddressing
{
    // Private data

        //- Reference to the octree
        const meshOctree& octree_;

        //- Reference to the dictionary
        const dictionary& meshDict_;

        //- Use DATA boxes
        bool useDATABoxes_;

        //- Number of created octree nodes
        mutable label nNodes_;

        //- Coordinates of octree nodes
        mutable pointField* octreePointsPtr_;

        //- Node labels
        mutable VRWGraph* nodeLabelsPtr_;

        //- Node leaves
        mutable FRWGraph<label, 8>* nodeLeavesPtr_;

        //- Identify which boxes should be used as mesh cells
        mutable List<direction>* boxTypePtr_;

        //- Identify created nodes as OUTERNODE or INNERNODE
        mutable List<direction>* nodeTypePtr_;


    // Additional addressing useful for some algorithms

        //- Faces of the octree
        mutable VRWGraph* octreeFacesPtr_;
        mutable labelLongList* octreeFacesOwnersPtr_;
        mutable labelLongList* octreeFacesNeighboursPtr_;

        //- Octree box-faces addressing
        mutable VRWGraph* leafFacesPtr_;

        //- Node-faces addressing
        mutable VRWGraph* nodeFacesPtr_;

        //- Leaf-leaves addressing
        mutable VRWGraph* leafLeavesPtr_;

        //- Edges of the octree
        mutable LongList<edge>* octreeEdgesPtr_;

        //- Edges-leaves addressing
        mutable VRWGraph* edgeLeavesPtr_;

        //- Leaf-edges addressing
        mutable VRWGraph* leafEdgesPtr_;

        //- Node-edges addressing
        mutable VRWGraph* nodeEdgesPtr_;

        //- Face-edges addressing
        mutable VRWGraph* faceEdgesPtr_;

        //- Edge-faces addressing
        mutable VRWGraph* edgeFacesPtr_;

        //- Global octree point label
        mutable labelLongList* globalPointLabelPtr_;

        //- Global point to local label addressing
        mutable Map<label>* globalPointToLocalPtr_;

        //- Point-processors addressing
        mutable VRWGraph* pointProcsPtr_;

        //- Global octree face label
        mutable labelLongList* globalFaceLabelPtr_;

        //- Global face label to local label addressing
        mutable Map<label>* globalFaceToLocalPtr_;

        //- Face-processors addressing
        mutable VRWGraph* faceProcsPtr_;

        //- Global leaf label
        mutable labelLongList* globalLeafLabelPtr_;

        //- Global leaf label to local label addressing for octree leaves
        mutable Map<label>* globalLeafToLocalPtr_;

        //- Leaf at procs
        mutable VRWGraph* leafAtProcsPtr_;


    // Private member functions

        //- Calculate octreePointsPtr_
        void createOctreePoints() const;

        //- Calculate nodeLabelsPtr_
        void createNodeLabels() const;

        //- Calculate nodeLeavesPtr_
        void createNodeLeaves() const;

        //- Assemble boxTypePtr_ list
        void findUsedBoxes() const;

        //- Calculate nodeTypePtr_
        void calculateNodeType() const;

        //- Calculate faces
        void createOctreeFaces() const;

        //- Calculate leaf-faces addressing
        void calculateLeafFaces() const;

        //- Calculate node-faces addressing
        void calculateNodeFaces() const;

        //- Calculate leaf-leaves addressing
        void calculateLeafLeaves() const;

        //- Calculate edges
        void createOctreeEdges() const;

        //- Calculate leaf-edges
        void calculateLeafEdges() const;

        //- Calculate edge-leaves addressing
        void calculateEdgeLeaves() const;

        //- Calculate edge-faces
        void calculateEdgeFaces() const;

        //- Clear allocated data
        void clearOut();
        void clearBoxTypes();
        void clearNodeAddressing();
        void clearOctreeFaces();
        void clearAddressing();
        void clearParallelAddressing();

        //- Check if distinct parts are glued together
        void checkGluedRegions();


    // Private member functions for parallel runs

        void calcGlobalPointLabels() const;
        void calcGlobalFaceLabels() const;
        void calcGlobalLeafLabels() const;


public:

    // Constructors

        //- Construct from surface and dictionary
        meshOctreeAddressing
        (
            const meshOctree& mo,
            const dictionary& dict,
            bool useDATABoxes = false
        );

        //- Disallow default bitwise copy construction
        meshOctreeAddressing(const meshOctreeAddressing&) = delete;


    //- Destructor
    ~meshOctreeAddressing();


    // Enumerators

        enum octreeCubeTypes
        {
            NONE = 0,
            MESHCELL = 1,
            BOUNDARY = 2,
            SPLITHEX = 4,
            INNERNODE = 8,
            OUTERNODE = 16,
            BOUNDARYNODE = 32
        };


    // Member Functions

        //- Check if there exist any non-signly connected edges and vertices
        //  in the octree mesh. This is primarily needed
        //  for the cartesian template
        void checkAndFixIrregularConnections();

        //- Return number of octree nodes
        inline label numberOfNodes() const;

        //- Return coordinates of octree vertices
        inline const pointField& octreePoints() const;

        //- Return nodeLabels
        inline const VRWGraph& nodeLabels() const;

        //- Return nodeLeaves
        inline const FRWGraph<label, 8>& nodeLeaves() const;

        //- Return which octree boxes are used for mesh creation
        inline const List<direction>& boxType() const;

        //- Return type of node (INNERNODE,or OUTERNODE)
        inline const List<direction>& nodeType() const;

        //- Set box type
        inline void setBoxType(const label boxI, const direction type);

        //- Return octree faces, created for MESHCELL boxes
        inline const VRWGraph& octreeFaces() const;

        //- Return owners of octree faces
        inline const labelLongList& octreeFaceOwner() const;

        //- Return neighbours of octree faces
        inline const labelLongList& octreeFaceNeighbour() const;

        //- Return octree box-faces addressing
        inline const VRWGraph& leafFaces() const;

        //- Checks if the face is intersected by the surface
        bool isIntersectedFace(const label fI) const;

        //- Return leaf-leaves addressing
        inline const VRWGraph& leafLeaves() const;

        //- Return node-faces addressing
        inline const VRWGraph& nodeFaces() const;

        //- Return octree edges, created for MESHCELL boxes
        inline const LongList<edge>& octreeEdges() const;

        //- Return edge-leaves addressing
        inline const VRWGraph& edgeLeaves() const;

        //- Return leaf-edges addressing
        inline const VRWGraph& leafEdges() const;

        //- Checks if the edge is intersected by the surface
        bool isIntersectedEdge(const label eI) const;
        void edgeIntersections
        (
            const label eI,
            DynList<point>& intersections
        ) const;

        //- Return node-edges addressing
        inline const VRWGraph& nodeEdges() const;

        //- Return face-edges addressing
        inline const VRWGraph& faceEdges() const;

        //- Return edge-faces addressing
        inline const VRWGraph& edgeFaces() const;

        //- Return const reference to meshOctree
        inline const meshOctree& octree() const;

        //- Find cubes around an edge (cubes must be at the same level)
        void cubesAroundEdge
        (
            const label leafI,
            const direction eI,
            FixedList<label, 4>& edgeCubes
        ) const;

        //- Find edge centre if it exists
        label findEdgeCentre(const label leafI, const direction eI) const;


    // Access to data needed for parallel execution

        //- Return global point labels
        inline const labelLongList& globalPointLabel() const;

        //- Global point label to local label. Only for processors points
        inline const Map<label>& globalToLocalPointAddressing() const;

        //- Processors which contain an octree point
        inline const VRWGraph& pointAtProcs() const;

        //- Return global labels of octree faces
        inline const labelLongList& globalFaceLabel() const;

        //- Return global face label to face label. Only for processor faces
        inline const Map<label>& globalToLocalFaceAddressing() const;

        //- Return processors which contain each octree face
        inline const VRWGraph& faceAtProcs() const;

        //- Return global labels of octree leaves
        inline const labelLongList& globalLeafLabel() const;

        //- Return processors which contain each octree leaf
        inline const VRWGraph& leafAtProcs() const;

        //- Return global leaf label to local label addressing
        inline const Map<label>& globalToLocalLeafAddressing() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshOctreeAddressing&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "utilities/octrees/meshOctree/meshOctreeAddressing/meshOctreeAddressingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
