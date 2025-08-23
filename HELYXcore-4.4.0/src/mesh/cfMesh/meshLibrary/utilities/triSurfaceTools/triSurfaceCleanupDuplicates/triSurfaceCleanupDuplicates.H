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
    triSurfaceCleanupDuplicates

Description
    Provides information regarding surface partitions on the surface
    triangulation, and connectivity between various surface partitions.

SourceFiles
    triSurfaceCleanupDuplicatesDuplicates.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceCleanupDuplicates_H
#define triSurfaceCleanupDuplicates_H

#include "utilities/meshes/triSurf/triSurf.H"
#include "utilities/containers/VRWGraph/VRWGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration
class meshOctree;

/*---------------------------------------------------------------------------*\
                Class triSurfaceCleanupDuplicates Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceCleanupDuplicates
{
    // Private data

        //- Tolerance for points
        const scalar tolerance_;

        //- Reference to triSurf
        triSurf& surf_;

        //- Reference to octree
        const meshOctree& octree_;

        //- New triangle labels in case some of them is removed
        labelLongList newTriangleLabel_;

        //- Check if the operation has already been performed
        bool done_;


    // Private member functions

        //- Check duplicate triangles
        bool checkDuplicateTriangles();

        //- Check duplicate points
        bool mergeDuplicatePoints();

        //- Update triangle mapping
        void updateTriangleLabels(const labelLongList&);


public:

    // Constructors

        //- Construct from meshOctree and tolerance
        triSurfaceCleanupDuplicates
        (
            const meshOctree& octree,
            const scalar tol=1e-6
        );

        //- Disallow default bitwise copy construction
        triSurfaceCleanupDuplicates
        (
            const triSurfaceCleanupDuplicates&
        ) = delete;


    //- Destructor
    ~triSurfaceCleanupDuplicates();


    // Member Functions

        //- Merge duplicate vertices and triangles
        void mergeIdentities();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceCleanupDuplicates&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
