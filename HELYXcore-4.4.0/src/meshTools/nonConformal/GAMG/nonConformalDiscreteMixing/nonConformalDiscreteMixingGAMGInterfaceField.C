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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "nonConformal/GAMG/nonConformalDiscreteMixing/nonConformalDiscreteMixingGAMGInterfaceField.H"
#include "nonConformal/polyPatches/nonConformalDiscreteMixing/intersection/discreteMixingIntersection.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalDiscreteMixingGAMGInterfaceField, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        nonConformalDiscreteMixingGAMGInterfaceField,
        lduInterface
    );
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        nonConformalDiscreteMixingGAMGInterfaceField,
        lduInterfaceField
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nonConformalDiscreteMixingGAMGInterfaceField::
nonConformalDiscreteMixingGAMGInterfaceField
(
    const GAMGInterface& GAMGNonConformalP,
    const lduInterfaceField& fineInterfaceField
)
:
    GAMGInterfaceField(GAMGNonConformalP, fineInterfaceField),
    nonConformalInterface_
    (
        refCast<const nonConformalDiscreteMixingGAMGInterface>
        (
            GAMGNonConformalP
        )
    )
{}


Foam::nonConformalDiscreteMixingGAMGInterfaceField::
nonConformalDiscreteMixingGAMGInterfaceField
(
    const GAMGInterface& GAMGNonConformalP,
    const int rank
)
:
    GAMGInterfaceField(GAMGNonConformalP, rank),
    nonConformalInterface_
    (
        refCast<const nonConformalDiscreteMixingGAMGInterface>
        (
            GAMGNonConformalP
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalDiscreteMixingGAMGInterfaceField::
~nonConformalDiscreteMixingGAMGInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nonConformalDiscreteMixingGAMGInterfaceField::updateInterfaceMatrix
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
        nonConformalInterface_.nbrPatch().interfaceInternalField(psiInternal)
    );

    if (nonConformalInterface_.owner())
    {
        pnf = nonConformalInterface_.intersection().interpolateToOwner(pnf);
    }
    else
    {
        pnf = nonConformalInterface_.nbrPatch().intersection()
                .interpolateToNeighbour(pnf);
    }

    this->addToInternalField(result, !add, coeffs, pnf);
}


// ************************************************************************* //
