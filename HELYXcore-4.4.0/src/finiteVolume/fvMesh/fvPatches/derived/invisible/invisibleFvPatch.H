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
    (c) 2011 OpenFOAM Foundation
    (c) 2020 Engys Ltd.

Class
    Foam::invisibleFvPatch

Description
    Foam::invisibleFvPatch

SourceFiles
    invisibleFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef invisibleFvPatch_H
#define invisibleFvPatch_H

#include "meshes/polyMesh/polyPatches/derived/invisible/invisiblePolyPatch.H"
#include "fvMesh/fvPatches/fvPatch/fvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class invisibleFvPatch Declaration
\*---------------------------------------------------------------------------*/

class invisibleFvPatch
:
    public fvPatch
{

public:

    //- Runtime type information
    TypeName(invisiblePolyPatch::typeName_());


    // Constructors

        //- Construct from components
        invisibleFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            fvPatch(patch, bm)
        {
            bool& act = this->active();
            act = false;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
