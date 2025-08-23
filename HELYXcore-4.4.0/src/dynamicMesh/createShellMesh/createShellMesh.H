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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::createShellMesh

Description
    Creates mesh by extruding a patch.

SourceFiles
    createShellMesh.C

    Extrudes into thickness direction.
    - bottom faces originate from reversed original faces (have turning index)
    - top faces originate from original faces (no turning index)

\*---------------------------------------------------------------------------*/

#ifndef createShellMesh_H
#define createShellMesh_H

#include "meshes/primitiveMesh/primitivePatch/primitiveFacePatch.H"
#include "containers/Lists/PackedList/PackedBoolList.H"
#include "primitives/Pair/labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyTopoChangeMap;
class polyTopoChange;
class globalMeshData;

/*---------------------------------------------------------------------------*\
                           Class createShellMesh Declaration
\*---------------------------------------------------------------------------*/

class createShellMesh
{
    // Private data

        //- Reference to patch to extrude
        const primitiveFacePatch& patch_;

        //- Per point on face the region
        const faceList& pointRegions_;

        //- From region to original patch point
        const labelList& regionPoints_;


        labelList cellToFaceMap_;

        labelList faceToFaceMap_;

        labelList faceToEdgeMap_;

        labelList pointToPointMap_;


    // Private Member Functions

        static void syncEdges
        (
            const globalMeshData&,
            const labelList&,
            const labelList&,
            const PackedBoolList& sameEdgeOrientation,
            const bool syncNonCollocated,

            PackedBoolList& isChangedEdge,
            DynamicList<label>& changedEdges,
            labelPairList& allEdgeData
        );


public:

    //- Runtime type information
    ClassName("createShellMesh");


    // Constructors

        //- Construct from mesh.
        createShellMesh
        (
            const primitiveFacePatch& patch,
            const faceList& pointRegions,
            const labelList& regionPoints
        );

        //- Disallow default bitwise copy construction
        createShellMesh(const createShellMesh&) = delete;


    // Member Functions


        // Access

            //- From region cell to patch face. Consecutively added so
            //  cell at layerI is at patchFacei*nLayers+layerI
            const labelList& cellToFaceMap() const
            {
                return cellToFaceMap_;
            }

            //- From region face to patch face. Contains turning index:
            //      > 0 : face in same orientation as patch face. face will
            //            be in top patch
            //      < 0 : face in opposite orientation as patch face. face will
            //            be in bottom patch
            //      = 0 : for all side faces
            const labelList& faceToFaceMap() const
            {
                return faceToFaceMap_;
            }

            //- From region side-face to patch edge. -1 for non-edge faces.
            const labelList& faceToEdgeMap() const
            {
                return faceToEdgeMap_;
            }

            //- From region point to patch point.
            const labelList& pointToPointMap() const
            {
                return pointToPointMap_;
            }


        // Other

            //- Helper: calculate point regions. The point region is the
            //  same on all faces connected to a point if they can be
            //  reached through a face-edge-face walk without crossing
            //  the nonManifoldEdge.
            //  syncNonCollocated = true: edges connected through cyclic become
            //  single region. false: kept separate.
            //  pointGlobalRegions : non-compact. Guaranteed to be the same
            //  across processors.
            //  pointLocalRegions : compact.
            //  localToGlobalRegion : local to global region.
            static void calcPointRegions
            (
                const globalMeshData& globalData,
                const primitiveFacePatch& patch,
                const PackedBoolList& nonManifoldEdge,
                const bool syncNonCollocated,
                faceList& pointGlobalRegions,
                faceList& pointLocalRegions,
                labelList& localToGlobalRegion
            );

            //- Play commands into polyTopoChange to create layer mesh.
            void setRefinement
            (
                const pointField& firstLayerThickness,
                const scalar expansionRatio,
                const label nLayers,
                const labelList& topPatchID,
                const labelList& bottomPatchID,
                const labelListList& extrudeEdgePatches,
                polyTopoChange& meshMod
            );

            //- Update any locally stored mesh information.
            void topoChange(const polyTopoChangeMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const createShellMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
