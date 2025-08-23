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

#include "offsetGAMGInterfaceField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "matrices/lduMatrix/lduMatrix/lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(offsetGAMGInterfaceField, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        offsetGAMGInterfaceField,
        lduInterface
    );
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        offsetGAMGInterfaceField,
        lduInterfaceField
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::offsetGAMGInterfaceField::offsetGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const lduInterfaceField& fineInterface
)
:
    GAMGInterfaceField(GAMGCp, fineInterface),
    slaveGAMGInterfaceField_
    (
        GAMGInterfaceField::New
        (
            refCast<const GAMGInterface>
            (
                refCast<const offsetLduInterfaceBase>(GAMGCp).slaveInterface()
            ),
            refCast<const offsetLduInterfaceFieldBase>
            (
                fineInterface
            ).slaveInterfaceField()
        )
    )
{}


Foam::offsetGAMGInterfaceField::offsetGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const int rank
)
:
    GAMGInterfaceField(GAMGCp, rank),
    slaveGAMGInterfaceField_
    (
        GAMGInterfaceField::New
        (
            refCast<const GAMGInterface>
            (
                refCast<const offsetLduInterfaceBase>(GAMGCp).slaveInterface()
            ),
            rank
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::offsetGAMGInterfaceField::~offsetGAMGInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::offsetGAMGInterfaceField::initInterfaceMatrixUpdate
(
    scalarField& result,
    const bool add,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    slaveGAMGInterfaceField_->initInterfaceMatrixUpdate
    (
        result,
        add,
        psiInternal,
        coeffs,
        cmpt,
        commsType
    );
}

void Foam::offsetGAMGInterfaceField::updateInterfaceMatrix
(
    scalarField& result,
    const bool add,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    slaveGAMGInterfaceField_->updateInterfaceMatrix
    (
        result,
        add,
        psiInternal,
        coeffs,
        cmpt,
        commsType
    );
}


// ************************************************************************* //
