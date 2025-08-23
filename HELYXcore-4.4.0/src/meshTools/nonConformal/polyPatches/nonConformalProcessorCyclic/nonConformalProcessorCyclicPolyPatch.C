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

\*---------------------------------------------------------------------------*/

#include "nonConformal/polyPatches/nonConformalProcessorCyclic/nonConformalProcessorCyclicPolyPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalProcessorCyclicPolyPatch, 0);
    addToRunTimeSelectionTable
    (
        polyPatch,
        nonConformalProcessorCyclicPolyPatch,
        dictionary
    );
}

// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::nonConformalProcessorCyclicPolyPatch::rename
(
    const wordList& newNames
)
{
    processorCyclicPolyPatch::rename(newNames);
    nonConformalCoupledPolyPatch::rename(newNames);
}


void Foam::nonConformalProcessorCyclicPolyPatch::reorder
(
    const labelUList& newToOldIndex
)
{
    processorCyclicPolyPatch::reorder(newToOldIndex);
    nonConformalCoupledPolyPatch::reorder(newToOldIndex);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nonConformalProcessorCyclicPolyPatch::
nonConformalProcessorCyclicPolyPatch
(
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const int myProcNo,
    const int neighbProcNo,
    const word& referPatchName,
    const word& origPatchName,
    const word& patchType
)
:
    processorCyclicPolyPatch
    (
        newName(referPatchName, myProcNo, neighbProcNo),
        size,
        start,
        index,
        bm,
        myProcNo,
        neighbProcNo,
        referPatchName,
        patchType
    ),
    nonConformalCoupledPolyPatch(*this, origPatchName)
{}


Foam::nonConformalProcessorCyclicPolyPatch::
nonConformalProcessorCyclicPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    processorCyclicPolyPatch(name, dict, index, bm, patchType),
    nonConformalCoupledPolyPatch(*this, dict)
{}


Foam::nonConformalProcessorCyclicPolyPatch::
nonConformalProcessorCyclicPolyPatch
(
    const nonConformalProcessorCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    processorCyclicPolyPatch(pp, bm),
    nonConformalCoupledPolyPatch(*this, pp)
{}


Foam::nonConformalProcessorCyclicPolyPatch::
nonConformalProcessorCyclicPolyPatch
(
    const nonConformalProcessorCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    processorCyclicPolyPatch(pp, bm, index, newSize, newStart),
    nonConformalCoupledPolyPatch(*this, pp)
{}


Foam::nonConformalProcessorCyclicPolyPatch::
nonConformalProcessorCyclicPolyPatch
(
    const nonConformalProcessorCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    processorCyclicPolyPatch(pp, bm, index, mapAddressing, newStart),
    nonConformalCoupledPolyPatch(*this, pp)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalProcessorCyclicPolyPatch::
~nonConformalProcessorCyclicPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nonConformalProcessorCyclicPolyPatch::write(Ostream& os) const
{
    processorCyclicPolyPatch::write(os);
    nonConformalCoupledPolyPatch::write(os);
}


// ************************************************************************* //
