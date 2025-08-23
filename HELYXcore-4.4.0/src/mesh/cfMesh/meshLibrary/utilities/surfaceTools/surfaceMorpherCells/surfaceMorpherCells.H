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
    surfaceMorpherCells

Description
    Changes the boundary of the mesh template, such that each boundary cell
    has only one boundary face.

SourceFiles
    surfaceMorpherCells.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceMorpherCells_H
#define surfaceMorpherCells_H

#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
#include "utilities/containers/Lists/labelLongList.H"
#include "primitives/bools/lists/boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class surfaceMorpherCells Declaration
\*---------------------------------------------------------------------------*/

class surfaceMorpherCells
{
    // Private data

        //- Mesh modifier
        polyMeshGen& mesh_;

        //- True for vertices which are at the boundary
        boolList boundaryVertex_;

        //- True for cells which are at the boundary
        List<direction> cellFlags_;

        //- New boundary faces and owner cells
        VRWGraph newBoundaryFaces_;
        labelLongList newBoundaryOwners_;
        labelLongList newBoundaryPatches_;


    // Private enumerators

        enum cellTypes
        {
            NONE = 0,
            BOUNDARY = 1,
            CHANGED = 2
        };


    // Private member functions

        //- Find surface vertices
        void findBoundaryVertices();

        //- Find boundary cells
        void findBoundaryCells();

        //- Remove cells with all vertices at the boundary
        bool removeCellsWithAllVerticesAtTheBoundary();

        //- Morph boundary faces
        bool morphBoundaryFaces();

        //- Morph internal faces
        bool morphInternalFaces();

        //- Replace the boundary of the mesh with the newly created one
        void replaceMeshBoundary();


public:

    // Constructors

        //- Construct from points and cell shapes
        surfaceMorpherCells(polyMeshGen& mesh);

        //- Disallow default bitwise copy construction
        surfaceMorpherCells(const surfaceMorpherCells&) = delete;


    //- Destructor
    ~surfaceMorpherCells();


    // Member Functions

        void morphMesh();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const surfaceMorpherCells&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
