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
    topologicalCleaner

Description
    Checks cells for invalid topologies. Decomposes bad cells into pyramids

SourceFiles
    topologicalCleaner.C

\*---------------------------------------------------------------------------*/

#ifndef topologicalCleaner_H
#define topologicalCleaner_H

#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
#include "primitives/bools/lists/boolList.H"
#include "containers/HashTables/HashSet/HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class topologicalCleaner Declaration
\*---------------------------------------------------------------------------*/

class topologicalCleaner
{
    // Private data

        polyMeshGen& mesh_;

        bool changed_;

        boolList decomposeCell_;


    // Private member functions

        //- Check and fix incorrectly connected faces
        void checkNonConsecutiveBoundaryVertices();

        //- Check and fix non-mappable faces
        void checkNonMappableFaces();

        //- Check and fix non-mappable cells
        void checkNonMappableCells();

        //- Decompose marked cells
        void decomposeCells();


public:

    // Constructors

        //- Construct from polyMeshGen
        topologicalCleaner(polyMeshGen& mesh);

        //- Disallow default bitwise copy construction
        topologicalCleaner(const topologicalCleaner&) = delete;


    //- Destructor
    ~topologicalCleaner();


    // Member Functions

        //- Clean topologically problematic cells
        bool cleanTopology();

        //- Check and fix vertices where two or more groups of cells meet
        void checkInvalidConnectionsForVerticesCells
        (
            labelHashSet* irregularNodes = nullptr
        );

        //- Check and fix vertices where two or more groups of faces meet
        void checkInvalidConnectionsForVerticesFaces
        (
            labelHashSet* irregularNodes = nullptr
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const topologicalCleaner&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
