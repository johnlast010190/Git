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
    meshOctreeModifier

Description
    Octree for mesh generation

SourceFiles
    meshOctreeModifier.C

\*---------------------------------------------------------------------------*/

#ifndef meshOctreeModifier_H
#define meshOctreeModifier_H

#include "utilities/octrees/meshOctree/meshOctree.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class meshOctreeModifier Declaration
\*---------------------------------------------------------------------------*/

class meshOctreeModifier
{
    // Private data

        //- Reference to the octree
        meshOctree& octree_;


    // Private member functions

        //- Correct refinement such that it produces 1-irregular octree
        void ensureCorrectRegularity(labelList& refineBox);

        //- Add additional refinement with the following rule. If a son
        //  of a father box is refined, then refine all other sons
        bool ensureCorrectRegularitySons(labelList& refineBox);


public:

    // Constructors

        //- Construct from octree
        meshOctreeModifier(meshOctree&);

        //- Disallow default bitwise copy construction
        meshOctreeModifier(const meshOctreeModifier&) = delete;


    //- Destructor
    ~meshOctreeModifier();


    // Member Functions

        //- Return octree
        inline const meshOctree& octree() const;

        //- Return neighbour processors
        inline labelList& neiProcsAccess();

        //- Return rootBox
        inline boundBox& rootBoxAccess();

        //- Return isRootInitialised_
        inline bool& isRootInitialisedAccess();

        //- Return search range
        inline scalar& searchRangeAccess();

        //- Return the range of leaves at each processor
        inline List<Pair<meshOctreeCubeCoordinates>>& neiRangeAccess();

        //- Return initial cube
        inline meshOctreeCube& initialCubeAccess();

        //- Return octree slots
        inline List<meshOctreeSlot>& dataSlotsAccess();

        //- Return the pointer to the meshOctreeCube at the given position
        inline meshOctreeCube* findCubeForPosition
        (
            const meshOctreeCubeCoordinates&
        ) const;

        //- Find leaves contained in a given boundBox
        inline void findLeavesContainedInBox
        (
            const boundBox&,
            DynList<const meshOctreeCube*, 256>&
        ) const;

        //- Return leaves
        inline LongList<meshOctreeCube*>& leavesAccess();

        //- Create leaves
        inline void createListOfLeaves();

        //- Mark additional layers around the leaves selected for refinement
        void markAdditionalLayers
        (
            labelList& refineBox,
            const label nLayers = 1
        ) const;

        //- Mark additional layers around the leaves selected for refinement
        void markAdditionalLayersOfFaceNeighbours
        (
            labelList& refineBox,
            const label nLayers = 1
        ) const;

        //- Mark additional layers around the leaves selected for refinement
        //  given on a box-by-box basis
        //  returns the number of boxes selected for refinement
        label markAdditionalLayers
        (
            labelList& refineBox,
            labelList& nLayers,
            List<direction>& targetRefLevel
        ) const;

        //- Refine leaves marked for refinement
        //  hexRefinement is activated when it is required to refine all
        //  sons of the same father, if a single son gets marked for refinement
        void refineSelectedBoxes
        (
            labelList& refineBox,
            const bool hexRefinement = false
        );

        //- Refine selected boxes and the boxes within the given range
        void refineSelectedBoxesAndAdditionalLayers
        (
            labelList& refineBox,
            const scalarList& refThickness
        );


    // Functions for parallel runs

        //- Distribute leaves of the initial octree to processors
        //  each processor creates a list of neighbouring processors
        void distributeLeavesToProcessors();

        //- Move octree cubes from one processor to another
        void loadDistribution(const direction usedType = 0);

        //- Refine the tree to add cubes transferred from other processors
        void refineTreeForCoordinates
        (
            const meshOctreeCubeCoordinates&,
            const short procNo = Pstream::myProcNo(),
            const direction cubeType = meshOctreeCubeBasic::UNKNOWN
        );

        void refineTreeForCoordinates
        (
            const meshOctreeCubeCoordinates&,
            const labelList& containedTriangles,
            const labelList& containedEdges,
            const short procNo = Pstream::myProcNo(),
            const direction cubeType = meshOctreeCubeBasic::UNKNOWN

        );

        //- Update the communication pattern between the processors
        //  this function is used for load balancing
        void updateCommunicationPattern();

        //- Add an additional layer of boxes which belongs to the neighbouring
        //  processors. This simplifies the process of extracting the mesh
        //  template from the octree in a parallel run
        void addLayerFromNeighbouringProcessors();

        //- Reduce memory consumption by deleting unnecessary data
        void reduceMemoryConsumption();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshOctreeModifier&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "utilities/octrees/meshOctree/meshOctreeModifier/meshOctreeModifierI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
