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

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Class
    polyMeshGenModifierAddCellByCell

Description
    Modifier for polyMeshGen

SourceFiles
    polyMeshGenModifierAddCellByCell.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshGenModifierAddCellByCell_H
#define polyMeshGenModifierAddCellByCell_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyMeshGenModifierAddCellByCell
: public polyMeshGenModifier
{
    // Private data
        //- number of faces
        label nFaces_;
        LongList<face> newFaces_;

        //- number of cells
        label nCells_;
        LongList<cell> newCells_;

public:

    // Constructors
        //- Construct from the reference to the mesh
        polyMeshGenModifierAddCellByCell(polyMeshGen& mesh);

    // Destructor
        ~polyMeshGenModifierAddCellByCell();

    // Member functions
        //- add a single cell (vertices must be added);
        void addCell(const faceList& cell);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
