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
    checkBoundaryFacesSharingTwoEdges

Description
    Checks cell connections at the surface. Cells with all vertices at
    the boundary are considered problematic because they get severly skewed
    after being mapped to the geometry surface in most cases. This code performs
    checks of cells with all vertices at the boundary with the rest of the mesh.
    Such cells having a neighbour with all internal faces are removed.

SourceFiles
    checkBoundaryFacesSharingTwoEdges.C

\*---------------------------------------------------------------------------*/

#ifndef checkBoundaryFacesSharingTwoEdges_H
#define checkBoundaryFacesSharingTwoEdges_H

#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
#include "containers/HashTables/HashSet/HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
            Class checkBoundaryFacesSharingTwoEdges Declaration
\*---------------------------------------------------------------------------*/

class checkBoundaryFacesSharingTwoEdges
{
    // Private data

        //- Reference to polyMeshGen
        polyMeshGen& mesh_;

        //- Pointer to meshSurfaceEngine
        mutable meshSurfaceEngine* meshSurfacePtr_;

        //- Number of boundary faces attached to a boundary vertex
        labelList nBndFacesAtBndPoint_;

        //- A list of boundary points which can be removed from the mesh
        boolList removeBndPoint_;


    // Private member functions

        //- Creates meshSurfaceEngine
        void createMeshSurface() const;

        //- Returns mesh surface
        inline const meshSurfaceEngine& meshSurface() const
        {
            if (!meshSurfacePtr_)
            {
                createMeshSurface();
            }

            return *meshSurfacePtr_;
        }

        //- Find the faces connected to each boundary edge
        //  mark the points which can be removed
        void findFacesAtBndEdge();

        //- Find the number of faces connected to the boundary vertex
        void findBndFacesAtBndVertex();

        //- Remove vertices from the faces
        void removeExcessiveVertices();

        //- Find boundary faces which shall be decomposed into triangles
        label findBndFacesForDecomposition(boolList& decomposeFace);


public:

    // Constructors

        //- Construct from polyMeshGen and meshOctree
        checkBoundaryFacesSharingTwoEdges(polyMeshGen& mesh);

        //- Disallow default bitwise copy construction
        checkBoundaryFacesSharingTwoEdges
        (
            const checkBoundaryFacesSharingTwoEdges&
        ) = delete;


    //- Destructor
    ~checkBoundaryFacesSharingTwoEdges();


    // Member Functions

        //- Find boundary points connected to two boundary faces, only
        void findPoints(labelHashSet& badPoints);

        //- Find boundary points connected to two boundary faces
        //  remove points if possible, decompose boundary faces otherwise
        bool improveTopology();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const checkBoundaryFacesSharingTwoEdges&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
