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
    Foam::boundaryCutter

Description
    Does modifications to boundary faces.

    Does
    - move boundary points
    - split boundary edges (multiple per edge if necessary)
    - face-centre decomposes boundary faces
    - diagonal split of boundary faces


SourceFiles
    boundaryCutter.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryCutter_H
#define boundaryCutter_H

#include "containers/HashTables/Map/Map.H"
#include "primitives/ints/lists/labelList.H"
#include "meshes/meshShapes/edge/EdgeMap.H"
#include "db/typeInfo/typeInfo.H"
#include "primitives/Pair/labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyTopoChange;
class polyTopoChangeMap;
class polyMesh;
class face;

/*---------------------------------------------------------------------------*\
                        Class boundaryCutter Declaration
\*---------------------------------------------------------------------------*/

class boundaryCutter
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Per edge sorted (start to end) list of points added.
        EdgeMap<labelList> edgeAddedPoints_;

        //- Per face the mid point added.
        Map<label> faceAddedPoint_;


    // Private Member Functions

        //- Get patch and zone info for face
        void getFaceInfo
        (
            const label facei,
            label& patchID,
            label& zoneID,
            label& zoneFlip
        ) const;

        //- Add cuts of edges to face
        face addEdgeCutsToFace(const label facei, const Map<labelList>&) const;

        //- Splits faces with multiple cut edges. Return true if anything split.
        bool splitFace
        (
            const label facei,
            const Map<point>& pointToPos,
            const Map<labelList>& edgeToAddedPoints,
            polyTopoChange& meshMod
        ) const;

        //- Add/modify facei for new vertices.
        void addFace
        (
            const label facei,
            const face& newFace,

            bool& modifiedFace,     // have we already 'used' facei?
            polyTopoChange& meshMod
        ) const;


public:

    //- Runtime type information
    ClassName("boundaryCutter");

    // Constructors

        //- Construct from mesh
        boundaryCutter(const polyMesh& mesh);


        //- Disallow default bitwise copy construction
        boundaryCutter(const boundaryCutter&) = delete;


    //- Destructor
    ~boundaryCutter();


    // Member Functions

        // Edit

            //- Do actual cutting with cut description. Inserts mesh changes
            //  into meshMod.
            //  pointToPos : new position for selected points
            //  edgeToCuts : per edge set of points that need to be introduced
            //  faceToSplit : per face the diagonal split
            //  faceToFeaturePoint : per face the feature point. Triangulation
            //                       around this feature point.
            void setRefinement
            (
                const Map<point>& pointToPos,
                const Map<List<point>>& edgeToCuts,
                const Map<labelPair>& faceToSplit,
                const Map<point>& faceToFeaturePoint,
                polyTopoChange& meshMod
            );

            //- Force recalculation of locally stored data on topological change
            void topoChange(const polyTopoChangeMap&);


        // Access

            //- Per edge a sorted list (start to end) of added points.
            const EdgeMap<labelList>& edgeAddedPoints() const
            {
                return edgeAddedPoints_;
            }

            //- Per face the mid point added.
            const Map<label>& faceAddedPoint() const
            {
                return faceAddedPoint_;
            }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const boundaryCutter&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
