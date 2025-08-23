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
    checkNonMappableCellConnections

Description
    Checks cell connections at the surface. Cells with all vertices at
    the boundary are considered problematic because they get severly skewed
    after being mapped to the geometry surface in most cases. This code performs
    checks of cells with all vertices at the boundary with the rest of the mesh.
    Such cells having a neighbour with all internal faces are removed.

SourceFiles
    checkNonMappableCellConnections.C

\*---------------------------------------------------------------------------*/

#ifndef checkNonMappableCellConnections_H
#define checkNonMappableCellConnections_H

#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
#include "containers/HashTables/HashSet/HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class checkNonMappableCellConnections Declaration
\*---------------------------------------------------------------------------*/

class checkNonMappableCellConnections
{
    // Private data

        //- Reference to polyMeshGen
        polyMeshGen& mesh_;

        //- type of cell
        labelList cellType_;


    // Enumerations

        enum cellTypes_
        {
            NONE = 0,
            INTERNALCELL = 1,
            BNDCELL = 2,
            ALLBNDVERTEXCELL = 4,
            INTERNALFACEGROUP = 8
        };


    // Private member functions

        //- Classify cells
        void findCellTypes();


public:

    // Constructors

        //- Construct from polyMeshGen and meshOctree
        checkNonMappableCellConnections(polyMeshGen& mesh);

        //- Disallow default bitwise copy construction
        checkNonMappableCellConnections
        (
            const checkNonMappableCellConnections&
        ) = delete;


    //- Destructor
    ~checkNonMappableCellConnections();


    // Member Functions

        //- Find problematic cells
        void findCells(labelHashSet& badCells);

        //- Find and remove problematic cells
        bool removeCells();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const checkNonMappableCellConnections&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
