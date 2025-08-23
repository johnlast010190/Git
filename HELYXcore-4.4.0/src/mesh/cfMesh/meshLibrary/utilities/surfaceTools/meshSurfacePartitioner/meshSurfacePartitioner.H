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
    meshSurfacePartitioner

Description
    Finds corners and edge points at the surface of the volume mesh

SourceFiles
    meshSurfacePartitioner.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfacePartitioner_H
#define meshSurfacePartitioner_H

#include "utilities/surfaceTools/meshSurfaceEngine/meshSurfaceEngine.H"
#include "utilities/containers/VRWGraph/VRWGraph.H"
#include "containers/HashTables/HashSet/HashSet.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class meshSurfacePartitioner Declaration
\*---------------------------------------------------------------------------*/

class meshSurfacePartitioner
{
    // Private data

        //- Reference to mesh surface engine
        const meshSurfaceEngine& meshSurface_;

        //- Pointer to face patches
        const labelList& facePatch_;

        //- Patches at a surface point
        VRWGraph pointPatches_;

        //- Labels of corner points
        labelHashSet corners_;

        //- Labels of edge points
        labelHashSet edgePoints_;

        //- Patch-patches addressing
        List<labelHashSet> patchPatches_;

        //- Number of edges attached to a surface point
        labelList nEdgesAtPoint_;

        //- Labels of feature edges
        labelHashSet featureEdges_;


    // Private member functions

        //- Find corners, edge points and addressing
        void calculateCornersEdgesAndAddressing();


public:

    // Constructors

        //- Construct from meshSurfaceEngine
        meshSurfacePartitioner(const meshSurfaceEngine&);

        //- Construct from meshSurfaceEngine and face patches
        meshSurfacePartitioner
        (
            const meshSurfaceEngine& meshSurface,
            const labelList& facePatch
        );

        //- Disallow default bitwise copy construction
        meshSurfacePartitioner(const meshSurfacePartitioner&) = delete;


    //- Destructor
    ~meshSurfacePartitioner();


    // Member Functions

        //- Return const reference to meshSurfaceEngine
        inline const meshSurfaceEngine& surfaceEngine() const
        {
            return meshSurface_;
        }

        //- Return const reference to a list containing boudary patch index
        //  for each boundary face
        inline const labelList& boundaryFacePatches() const
        {
            return facePatch_;
        }

        //- Return const reference to a graph of patches attached to a surface
        //  point
        inline const VRWGraph& pointPatches() const
        {
            return pointPatches_;
        }

        //- Return labels of corner points (from the list of boundary points)
        inline const labelHashSet& corners() const
        {
            return corners_;
        }

        //- Return labels of edge points (from the list of boundary points)
        inline const labelHashSet& edgePoints() const
        {
            return edgePoints_;
        }

        //- Return the number of feature edges attached to a boundary point
        inline label numberOfFeatureEdgesAtPoint(const label bpI) const
        {
            return nEdgesAtPoint_[bpI];
        }

        //- Return patch-patches addressing
        inline const List<labelHashSet>& patchPatches() const
        {
            return patchPatches_;
        }

        //- return labels of boundary edges which are feature edges
        inline const labelHashSet& featureEdges() const
        {
            return featureEdges_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshSurfacePartitioner&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
