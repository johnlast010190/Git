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
    triSurfacePartitioner

Description
    Provides information regarding surface partitions on the surface
    triangulation, and connectivity between various surface partitions.

SourceFiles
    triSurfacePartitioner.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfacePartitioner_H
#define triSurfacePartitioner_H

#include "utilities/meshes/triSurf/triSurf.H"
#include "utilities/containers/DynList/DynList.H"
#include "containers/HashTables/HashSet/HashSet.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class triSurfacePartitioner Declaration
\*---------------------------------------------------------------------------*/

class triSurfacePartitioner
{
    // Private data

        //- Reference to triSurf
        const triSurf& surface_;

        //- Corner nodes
        labelList corners_;
        List<DynList<label>> cornerPatches_;

        //- Information which partitions share an edge with a given partition
        List<labelHashSet> patchPatches_;

        //- Edge partitions
        labelList edgeGroups_;

        //- Information which edge groups share a corner
        //  with a given group
        List<labelHashSet> edgeGroupEdgeGroups_;

        //- Edge groups between surface patches
        std::map<std::pair<label, label>, labelHashSet> patchesEdgeGroups_;

        //- Corners shared by edge groups
        std::map<std::pair<label, label>, labelHashSet> edgeGroupsCorners_;


    // Private member functions

        //- Calculate patch addressing
        void calculatePatchAddressing();

        //- Find surface corners
        void calculateCornersAndAddressing();

        //- Calculate patch-patches addressing
        void calculatePatchPatches();

        //- Calculate edge groups
        void calculateEdgeGroups();

        //- Calculate surface patch to edge groups addressing
        void calculatePatchToEdgeGroups();

        //- Calculate edge groups to corner addressing
        void calculateEdgeGroupsToCorners();


public:

    // Constructors

        //- Construct from triSurf
        triSurfacePartitioner(const triSurf& surface);

        //- Disallow default bitwise copy construction
        triSurfacePartitioner(const triSurfacePartitioner&) = delete;


    //- Destructor
    ~triSurfacePartitioner();


    // Member Functions

        //- Return corner nodes
        const labelList& corners() const;

        //- Return corner patches
        const List<DynList<label>>& cornerPatches() const;

        //- Return patch-patches addressing
        const List<labelHashSet>& patchPatches() const;

        //- Return edge groups. Edges which are not
        //  feature edges are set to -1
        const labelList& edgeGroups() const;

        //- Edge group - edge groups addressing
        const List<labelHashSet>& edgeGroupEdgeGroups() const;

        //- Return groups of feature edges shared by the given patches
        void edgeGroupsSharedByPatches
        (
            const label partition1,
            const label partition2,
            DynList<label>& edgePartitions
        ) const;

        //- Return corners shared shared by the given edge groups
        void cornersSharedByEdgeGroups
        (
            const label edgePartition1,
            const label edgePartition2,
            DynList<label>& corners
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const triSurfacePartitioner&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
