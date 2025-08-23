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
    meshSurfaceEdgeExtractor

Description
    Extracts edges and corners to cornform the mesh to the geometry

SourceFiles
    meshSurfaceEdgeExtractor.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceEdgeExtractor_H
#define meshSurfaceEdgeExtractor_H

#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
#include "primitives/bools/lists/boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
                Class meshSurfaceEdgeExtractor Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceEdgeExtractor
{
    // Private data

        //- Mesh
        polyMeshGen& mesh_;

        label nPoints_;

        boolList boundaryCell_;

        List<direction> nFacesInCell_;

        //- Reference to the octree
        const meshOctree& meshOctree_;

        //- Regions for boundary vertices
        VRWGraph pointRegions_;


    // Private member functions

        //- Remove old boundary faces
        void removeOldBoundaryFaces();

        //- Create vertices on surface edges
        void createEdgeVertices();

        //- Create new boundary faces
        void createBoundaryFaces();


public:

    // Constructors

        //- Construct from mesh data
        meshSurfaceEdgeExtractor
        (
            polyMeshGen& mesh,
            const meshOctree& octree,
            const labelList& pointRegion
        );

        //- Disallow default bitwise copy construction
        meshSurfaceEdgeExtractor(const meshSurfaceEdgeExtractor&) = delete;


    // Destructor
    ~meshSurfaceEdgeExtractor();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceEdgeExtractor&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
