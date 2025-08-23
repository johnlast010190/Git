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
    Foam::nonConformalProcessorCyclicPointPatch

Description
    Constraint patch for nonConformalProcessorCyclic couplings.

SourceFiles
    nonConformalProcessorCyclicPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalProcessorCyclicPointPatch_H
#define nonConformalProcessorCyclicPointPatch_H

#include "meshes/pointMesh/pointPatches/constraint/processorCyclic/processorCyclicPointPatch.H"
#include "nonConformal/polyPatches/nonConformalProcessorCyclic/nonConformalProcessorCyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class nonConformalProcessorCyclicPointPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalProcessorCyclicPointPatch
:
    public processorCyclicPointPatch
{
public:

    //- Runtime type information
    TypeName(nonConformalProcessorCyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        nonConformalProcessorCyclicPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        )
        :
            processorCyclicPointPatch(patch, bm)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
