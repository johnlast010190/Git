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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "matrices/lduMatrix/solvers/GAMG/interfaceFields/cyclicGAMGInterfaceField/cyclicGAMGInterfaceField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "matrices/lduMatrix/lduMatrix/lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicGAMGInterfaceField, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        cyclicGAMGInterfaceField,
        lduInterface
    );
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        cyclicGAMGInterfaceField,
        lduInterfaceField
    );

    // Add under name cyclicSlip
    addNamedToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        cyclicGAMGInterfaceField,
        lduInterface,
        cyclicSlip
    );
    addNamedToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        cyclicGAMGInterfaceField,
        lduInterfaceField,
        cyclicSlip
    );

    // Add under name nonConformalCyclic
    addNamedToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        cyclicGAMGInterfaceField,
        lduInterface,
        nonConformalCyclic
    );
    addNamedToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        cyclicGAMGInterfaceField,
        lduInterfaceField,
        nonConformalCyclic
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cyclicGAMGInterfaceField::cyclicGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const lduInterfaceField& fineInterface
)
:
    GAMGInterfaceField(GAMGCp, fineInterface),
    cyclicInterface_(refCast<const cyclicGAMGInterface>(GAMGCp)),
    rank_(0)
{
    const cyclicLduInterfaceField& p =
        refCast<const cyclicLduInterfaceField>(fineInterface);

    rank_ = p.rank();
}


Foam::cyclicGAMGInterfaceField::cyclicGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const int rank
)
:
    GAMGInterfaceField(GAMGCp, rank),
    cyclicInterface_(refCast<const cyclicGAMGInterface>(GAMGCp)),
    rank_(rank)
{}


// * * * * * * * * * * * * * * * * Desstructor * * * * * * * * * * * * * * * //

Foam::cyclicGAMGInterfaceField::~cyclicGAMGInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cyclicGAMGInterfaceField::updateInterfaceMatrix
(
    scalarField& result,
    const bool add,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    // Get neighbouring field
    scalarField pnf
    (
        cyclicInterface_.nbrPatch().interfaceInternalField(psiInternal)
    );

    transformCoupleField(pnf, cmpt);

    this->addToInternalField(result, !add, coeffs, pnf);
}


// ************************************************************************* //
