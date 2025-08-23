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
    renameBoundaryPatches

Description
    Renames boundary patches based on user intput and sets correct patch types

SourceFiles
    renameBoundaryPatches.C

\*---------------------------------------------------------------------------*/

#ifndef renameBoundaryPatches_H
#define renameBoundaryPatches_H

#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class IOdictionary;

/*---------------------------------------------------------------------------*\
                           Class renameBoundaryPatches Declaration
\*---------------------------------------------------------------------------*/

class renameBoundaryPatches
{
    // Private data

        //- Mesh
        polyMeshGen& mesh_;

        //- Dictionary containing relevant information
        const IOdictionary& meshDict_;


    // Private member functions

        //- Calculate new boundary and replace the existing one
        void calculateNewBoundary();


public:

    // Constructors

        //- Construct from mesh and IOdictionary
        renameBoundaryPatches
        (
            polyMeshGen& mesh,
            const IOdictionary& meshDict
        );

        //- Disallow default bitwise copy construction
        renameBoundaryPatches(const renameBoundaryPatches&) = delete;


    //- Destructor
    ~renameBoundaryPatches();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const renameBoundaryPatches&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
