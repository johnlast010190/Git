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
    Foam::removePoints

Description
    Removes selected points from mesh and updates faces using these
    points.

SourceFiles
    removePoints.C

\*---------------------------------------------------------------------------*/

#ifndef removePoints_H
#define removePoints_H

#include "db/typeInfo/typeInfo.H"
#include "primitives/bools/lists/boolList.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "meshes/meshShapes/face/faceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyTopoChange;
class polyTopoChangeMap;
class face;

/*---------------------------------------------------------------------------*\
                           Class removePoints Declaration
\*---------------------------------------------------------------------------*/

class removePoints
{

    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Whether undoable
        const bool undoable_;

        //- If undoable: deleted points
        pointField savedPoints_;

        //- If undoable: per stored face the original mesh face label
        labelList savedFaceLabels_;

        //- If undoable: per stored face the vertices. Negative indices
        //  refer to deletedPoints_
        faceList savedFaces_;


    // Private Member Functions

        //- Change the vertices of the face whilst keeping everything else
        //  (patch, zone) the same.
        void modifyFace
        (
            const label facei,
            const face&,
            polyTopoChange&
        ) const;


public:

    //- Runtime type information
    ClassName("removePoints");


    // Constructors

        //- Construct from mesh
        removePoints(const polyMesh& mesh, const bool undoable = false);

        //- Disallow default bitwise copy construction
        removePoints(const removePoints&) = delete;


    // Member Functions

            //- If undoable: affected face labels. Already restored faces
            //  will be -1.
            const labelList& savedFaceLabels() const
            {
                return savedFaceLabels_;
            }


        // Helper functions

            //- Mark in pointCanBeDeleted the points that can be deleted
            //  (parallel synchronised) and returns the global number of these
            //  points. (this number is the global number before synchronisation
            //  so might be off!)
            //  A point can be deleted if
            //  - it is not used by any edge.
            //  or
            //  - is used by only two edges
            //  - these two edges are sufficiently in line (cos > minCos)
            //  - all processors agree that point can be deleted.
            label countPointUsage
            (
                const scalar minCos,
                boolList& pointCanBeDeleted
            ) const;

        // Topology changes

            //- Play commands into polyTopoChange to remove points. Gets
            //  boolList (output of countPointUsage) as input.
            //  Does no check for whether resulting face is legal.
            //  Since pointCanBeDeleted is synced all coupled faces should
            //  decide the same.
            void setRefinement(const boolList&, polyTopoChange&);

            //- Force recalculation of locally stored data on topological change
            void topoChange(const polyTopoChangeMap&);

            //- Given set of faces to restore calculates a consistent set of
            //  saved faces (indices into savedFaces_) and saved vertices
            //  (indices into savedPoints_) to restore. The undoFaces have to
            //  be synced.
            void getUnrefimentSet
            (
                const labelList& undoFaces,
                labelList& localFaces,
                labelList& localPoints
            ) const;

            //- Restore selected faces and vertices.
            void setUnrefinement
            (
                const labelList& localFaces,
                const labelList& localPoints,
                polyTopoChange&
            );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const removePoints&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
