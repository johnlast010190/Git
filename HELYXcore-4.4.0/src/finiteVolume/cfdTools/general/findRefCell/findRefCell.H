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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2024 Engys Ltd

InNamespace
    Foam

Description
    Find the reference cell nearest (in index) to the given cell but which
    is not on a cyclic, symmetry or processor patch.

SourceFiles
    findRefCell.C

\*---------------------------------------------------------------------------*/

#ifndef findRefCell_H
#define findRefCell_H

#include "primitives/ints/label/label.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Find first processor which has a nonzero number of cells and return index
//  zero on that processor and -1 on the others. Returns (cellNo, processorNo).
Pair<label> getFirstCell(const fvMesh& mesh);


//- If the field fieldRef needs referencing find the reference cell nearest
//  (in index) to the given cell looked-up for field, but which is not on a
//  cyclic, symmetry or processor patch and return true, otherwise return false.
bool setRefCell
(
    const volScalarField& field,
    const volScalarField& fieldRef,
    const dictionary& dict,
    label& refCelli,
    scalar& refValue,
    const bool forceReference = false
);

//- If the field needs referencing find the reference cell nearest
//  (in index) to the given cell looked-up for field, but which is not on a
//  cyclic, symmetry or processor patch and return true, otherwise return false.
bool setRefCell
(
    const volScalarField& field,
    const dictionary& dict,
    label& refCelli,
    scalar& refValue,
    const bool forceReference = false
);

//- Return the current value of field in the reference cell
scalar getRefCellValue
(
    const volScalarField& field,
    const label refCelli
);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
