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
    voronoiMeshExtractor

Description
    A class which extracts tet mesh out of an octree structure

SourceFiles
    voronoiMeshExtractor.C

\*---------------------------------------------------------------------------*/

#ifndef voronoiMeshExtractor_H
#define voronoiMeshExtractor_H

#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
#include "utilities/meshes/primitives/partTet/partTet.H"
#include "utilities/tetrahedra/tetCreatorOctree/tetCreatorOctree.H"
#include "utilities/containers/VRWGraph/VRWGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class voronoiMeshExtractor Declaration
\*---------------------------------------------------------------------------*/

class voronoiMeshExtractor
{
    // Private data

        //- Create tets
        tetCreatorOctree tetCreator_;

        //- Reference to the mesh
        polyMeshGen& mesh_;

        //- Edges connected to a vertex in the tet mesh
        mutable VRWGraph* pointEdgesPtr_;

        //- Edges of the tetrahedrization
        mutable LongList<edge>* edgesPtr_;

        //- Tets around each edge
        mutable VRWGraph* edgeTetsPtr_;

        //- Is edge at the boundary
        mutable boolList* boundaryEdgePtr_;


    // Private member functions

        //- Create and return addressing
        void createAddressing() const;

        const VRWGraph& pointEdges() const;
        const LongList<edge>& edges() const;
        const VRWGraph& edgeTets() const;
        const boolList& boundaryEdge() const;

        //- Create points of the voronoi mesh
        void createPoints();

        //- Create mesh data
        void createPolyMesh();

        //- Delete all pointer data
        void clearOut();


    // Private static data

        //- Positive orientation
        static label sameOrientation_[6];

        //- Opposite orientation
        static label oppositeOrientation_[6];


public:

    // Constructors

        //- Construct from octree and mesh data
        voronoiMeshExtractor
        (
            const meshOctree& octree,
            const IOdictionary& meshDict,
            polyMeshGen& mesh
        );

        //- Disallow default bitwise copy construction
        voronoiMeshExtractor(const voronoiMeshExtractor&) = delete;


    //- Destructor
    ~voronoiMeshExtractor();


    // Member Functions

        //- Create mesh
        void createMesh();


    // Member Operators


        //- Disallow default bitwise assignment
        void operator=(const voronoiMeshExtractor&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
