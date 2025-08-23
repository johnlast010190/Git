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
    meshSurfaceEdgeExtractorFUN

Description
    Stores boundary faces into patches and captures edges and corners
    by inserting fundamental mesh sheets

SourceFiles
    meshSurfaceEdgeExtractorFUN.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceEdgeExtractorFUN_H
#define meshSurfaceEdgeExtractorFUN_H

#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
#include "primitives/bools/lists/boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
                Class meshSurfaceEdgeExtractorFUN Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceEdgeExtractorFUN
{
    // Private data

        //- Mesh
        polyMeshGen& mesh_;

        //- Octree
        const meshOctree& meshOctree_;

        //- Mesh surface pointer
        meshSurfaceEngine* surfaceEnginePtr_;

        //- Shall the procedure generate an initial wrapper sheet
        const bool createWrapperSheet_;


    // Private member functions

        //- Return reference to surface engine
        meshSurfaceEngine& surfaceEngine();

        //- Clear mesh suirface engine
        void clearOut();

        //- Distribute boundary faces into patches
        void distributeBoundaryFaces();

        //- Check whether all corners in the surface mesh are present
        //  in the volume mesh
        void reviseCorners();

        //- Check whether all edges in the surface mesh are present
        //  in the volume mesh
        void reviseEdges();

        //- Find corner vertices and correct patches
        void findCornersAndCorrectPatches();

        //- Create fundamental sheets
        void createBasicFundamentalSheets();

        //- Smooth the surface
        void smoothMeshSurface();

        //- Modify fundamental sheets
        void improveQualityOfFundamentalSheets();

        //- Re-map points after edges have been extracted
        void remapBoundaryPoints();


public:

    // Constructors

        //- Construct from mesh data
        meshSurfaceEdgeExtractorFUN
        (
            polyMeshGen& mesh,
            const meshOctree& octree,
            const bool createWrapperSheet = true
        );

        //- Disallow default construct
        meshSurfaceEdgeExtractorFUN() = delete;

        //- Disallow default bitwise copy construction
        meshSurfaceEdgeExtractorFUN
        (
            const meshSurfaceEdgeExtractorFUN&
        ) = delete;


    //- Destructor
    ~meshSurfaceEdgeExtractorFUN();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceEdgeExtractorFUN&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
