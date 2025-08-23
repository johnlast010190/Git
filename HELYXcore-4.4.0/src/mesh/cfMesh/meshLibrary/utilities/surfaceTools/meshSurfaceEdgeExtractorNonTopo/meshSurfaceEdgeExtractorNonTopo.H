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

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Class
    meshSurfaceEdgeExtractorNonTopo

Description
    Stores boundary faces into patches and captures edges and corners

SourceFiles
    meshSurfaceEdgeExtractorNonTopo.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceEdgeExtractorNonTopo_H
#define meshSurfaceEdgeExtractorNonTopo_H

#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;

/*---------------------------------------------------------------------------*\
                Class meshSurfaceEdgeExtractorNonTopo Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceEdgeExtractorNonTopo
{
    // Private data

        //- Reference to polyMeshGen
        polyMeshGen& mesh_;

        //- Reference to meshOctree
        const meshOctree& meshOctree_;


    // Private member functions

        //- Decompose faces with more than one edge at the same feature edge
        void decomposeBoundaryFaces();

        //- Re-map points after edges have been extracted
        void remapBoundaryPoints();


public:

    // Constructors

        //- Construct from mesh data
        meshSurfaceEdgeExtractorNonTopo
        (
            polyMeshGen& mesh,
            const meshOctree& octree
        );

        //- Disallow default bitwise copy construction
        meshSurfaceEdgeExtractorNonTopo
        (
            const meshSurfaceEdgeExtractorNonTopo&
        ) = delete;


    //- Destructor
    ~meshSurfaceEdgeExtractorNonTopo();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceEdgeExtractorNonTopo&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
