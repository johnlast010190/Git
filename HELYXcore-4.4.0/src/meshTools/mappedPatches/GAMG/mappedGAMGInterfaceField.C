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
    (c) 2010-2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "mappedPatches/GAMG/mappedGAMGInterfaceField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "matrices/lduMatrix/lduMatrix/lduMatrix.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedGAMGInterfaceField, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        mappedGAMGInterfaceField,
        lduInterface
    );
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        mappedGAMGInterfaceField,
        lduInterfaceField
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mappedGAMGInterfaceField::mappedGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const lduInterfaceField& fineInterface
)
:
    GAMGInterfaceField(GAMGCp, fineInterface),
    regionCoupledGAMGInterface_
    (
        refCast<const mappedGAMGInterface>(GAMGCp)
    )
{}


Foam::mappedGAMGInterfaceField::mappedGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const int rank
)
:
    GAMGInterfaceField(GAMGCp, rank),
    regionCoupledGAMGInterface_
    (
        refCast<const mappedGAMGInterface>(GAMGCp)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedGAMGInterfaceField::~mappedGAMGInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mappedGAMGInterfaceField::updateInterfaceMatrix
(
    scalarField& result,
    const bool add,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    if (regionCoupledGAMGInterface_.nbrPatchID() == -1)
    {
        return;
    }

    // Get neighbouring field
    scalarField pnf
    (
        regionCoupledGAMGInterface_.nbrPatch().interfaceInternalField(psiInternal)
    );

    // Transform according to the transformation tensors
    //transformCoupleField(pnf, cmpt); //TODO: Need to allow for transformations

    regionCoupledGAMGInterface_.interpolate(pnf);

    const labelUList& faceCells = regionCoupledGAMGInterface_.faceCells();

    if (add)
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
        }
    }
    else
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] += coeffs[elemI]*pnf[elemI];
        }
    }
}



// ************************************************************************* //
