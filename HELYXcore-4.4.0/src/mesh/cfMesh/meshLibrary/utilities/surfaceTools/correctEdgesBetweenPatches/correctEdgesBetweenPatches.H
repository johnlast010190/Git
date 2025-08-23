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
    correctEdgesBetweenPatches

Description
    Changes mesh surface until every boundary face shares only one edge with
    the faces in other boundary patch

SourceFiles
    correctEdgesBetweenPatches.C

\*---------------------------------------------------------------------------*/

#ifndef correctEdgesBetweenPatches_H
#define correctEdgesBetweenPatches_H

#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
#include "utilities/containers/Lists/labelLongList.H"
#include "primitives/bools/lists/boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
                Class correctEdgesBetweenPatches Declaration
\*---------------------------------------------------------------------------*/

class correctEdgesBetweenPatches
{
    // Private data

        //- Mesh
        polyMeshGen& mesh_;

        //- Pointer to surface engine
        mutable const meshSurfaceEngine* msePtr_;

        //- Boundary data
        wordList patchNames_;
        wordList patchTypes_;
        VRWGraph newBoundaryFaces_;
        labelLongList newBoundaryOwners_;
        labelLongList newBoundaryPatches_;

        //- Holds data which cells have to be decomposed
        boolList decomposeCell_;
        bool decompose_;


    // Private member functions

        //- Construct and return mesh surface
        const meshSurfaceEngine& meshSurface() const;

        //- Delete mesh surface
        void clearMeshSurface();

         //- Replace boundary
        void replaceBoundary();

        //- Decompose corrected cells
        void decomposeCorrectedCells();

        //- Decompose bnd faces for which the feature edges are not connected
        //  over vertices and internal faces which have more than one feature
        //  edge at the boundary
        void decomposeProblematicFaces();

        //- Decompose cells with faces at concave boundary edges
        void decomposeConcaveFaces();

        //- Perform decompose remaining faces having more than one feature edge
        void patchCorrection();


public:

    // Constructors

        //- Construct from mesh data
        correctEdgesBetweenPatches(polyMeshGen& mesh);

        //- Disallow default bitwise copy construction
        correctEdgesBetweenPatches(const correctEdgesBetweenPatches&) = delete;


    //- Destructor
    ~correctEdgesBetweenPatches();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const correctEdgesBetweenPatches&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
