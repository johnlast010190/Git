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
    tetCreatorOctree

Description
    A class which creates tets from the octree structure

SourceFiles
    tetCreatorOctree.C

\*---------------------------------------------------------------------------*/

#ifndef tetCreatorOctree_H
#define tetCreatorOctree_H

#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
#include "utilities/meshes/primitives/partTet/partTet.H"
#include "utilities/octrees/meshOctree/meshOctreeAddressing/meshOctreeAddressing.H"
#include "utilities/octrees/meshOctree/meshOctreeCube/meshOctreeCubeCoordinates.H"
#include "utilities/containers/VRWGraph/VRWGraph.H"
#include "utilities/containers/Lists/labelLongList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class tetCreatorOctree Declaration
\*---------------------------------------------------------------------------*/

class tetCreatorOctree
{
    // Private data

        //- Helper for searching coordinates of cubes around an edge
        static const meshOctreeCubeCoordinates edgeCoordinates_[12][4];

        //- Helper for finding face centres of cubes sharing an edge
        static const label faceCentreHelper_[3][4];

        //- Reference to the octree
        meshOctreeAddressing octreeCheck_;

        //- Points of the tetrahedrisation
        LongList<point> tetPoints_;

        //- Tetrahedra
        LongList<partTet> tets_;

        //- Octree leaves sorted according to their level
        List<labelLongList> sortedLeaves_;

        //- Node labels of vertices created inside split-hex boxes
        VRWGraph* subNodeLabelsPtr_;

        //- Cube centre label
        labelList* cubeLabelPtr_;

        //- Cube face label
        VRWGraph* faceCentreLabelPtr_;

        //- Are tets created or not
        bool created_;


    // Private member functions

        //- Deletes all pointer data
        void clearOut();

        //- Find elements which will be used as mesh cells
        void selectElements();


    // Private member functions for creating tetrahedra

        //- Create tetPoints_ and necessary addressing
        void createPointsAndAddressing();

        //- Create faceCentreLabelPtr_
        void createFaceCentreLabels();

        //- Create tetrahedra from faces, owner and neighbour
        void createTetsAroundEdges();

        //- Create tetrahedra from split faces
        void createTetsFromSplitFaces();

        //- Create tetrahedra from faces which share at least one edge
        //  with a refined cube
        void createTetsFromFacesWithCentreNode();

        //- Create tetrahedra from faces with split edges
        void createTetsAroundSplitEdges();

        //- Helper funtion which checks validity of a created tet and appends
        //  it if valid
        inline void checkAndAppendTet(const partTet);

        //- Function containing the workflow
        void createTets();


public:

    // Constructors

        //- Construct from octree and IOdictionary
        tetCreatorOctree
        (
            const meshOctree& octree,
            const IOdictionary& meshDict
        );

        //- Disallow default bitwise copy construction
        tetCreatorOctree(const tetCreatorOctree&) = delete;


    //- Destructor
    ~tetCreatorOctree();


    // Member Functions

        const LongList<point>& tetPoints() const
        {
            if (!created_)
                FatalErrorInFunction
                    << "Tets are not created!"
                    << exit(FatalError);

            return tetPoints_;
        }

        const LongList<partTet>& tets() const
        {
            if (!created_)
                FatalErrorInFunction
                    << "Tets are not created!"
                    << exit(FatalError);

            return tets_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const tetCreatorOctree&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "utilities/tetrahedra/tetCreatorOctree/tetCreatorOctreeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
