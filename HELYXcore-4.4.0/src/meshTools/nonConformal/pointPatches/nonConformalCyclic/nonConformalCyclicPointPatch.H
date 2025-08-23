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
    (c) 2021-2022 OpenFOAM Foundation
    (c) 2022-2025 Engys Ltd.

Class
    Foam::nonConformalCyclicPointPatch

Description
    Constraint patch for nonConformalCyclic couplings.

SourceFiles
    nonConformalCyclicPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalCyclicPointPatch_H
#define nonConformalCyclicPointPatch_H

#include "meshes/pointMesh/pointPatches/constraint/cyclic/cyclicPointPatch.H"
#include "nonConformal/polyPatches/nonConformalCyclic/nonConformalCyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class nonConformalCyclicPointPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalCyclicPointPatch
:
    public cyclicPointPatch
{
public:

    //- Runtime type information
    TypeName(nonConformalCyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        nonConformalCyclicPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        )
        :
            cyclicPointPatch(patch, bm)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
