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
    (c) 2013-2016 OpenFOAM Foundation
    (c) 2020-2024 Engys Ltd.

Class
    Foam::orientFaceZone

Description
    Corrects orientation of faceZone.

    - correct in parallel - excludes coupled faceZones from walk
    - correct for non-manifold faceZones - restarts walk

SourceFiles
    orientFaceZone.C

\*---------------------------------------------------------------------------*/

#ifndef orientFaceZone_H
#define orientFaceZone_H

#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class orientFaceZone Declaration
\*---------------------------------------------------------------------------*/

class orientFaceZone
{
    // Private data

        //- polyMesh reference
        polyMesh& mesh_;

        //- faceZone name
        const word& fzName_;


public:

    //- Runtime type information
    ClassName("orientFaceZone");


    // Constructors

        //- Construct from cells to cut and pattern of cuts
        orientFaceZone
        (
            polyMesh&,
            const word&
        );

        //- Disallow default bitwise copy construction
        orientFaceZone(const orientFaceZone&) = delete;


    //- Destructor
    ~orientFaceZone(){};


    // Member Functions

        //- Create and write mesh surface
        void orientByPoint(const point&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const orientFaceZone&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
