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
    (c) 2012-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "derivedFvPatchFields/energyJump/energyJumpAMI/energyJumpAMIFvPatchScalarField.H"
#include "fields/fvPatchFields/derived/fixedJumpAMI/fixedJumpAMIFvPatchFields.H"
#include "basicThermo/basicThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::energyJumpAMIFvPatchScalarField::energyJumpAMIFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpAMIFvPatchField<scalar>(p, iF)
{}


Foam::energyJumpAMIFvPatchScalarField::energyJumpAMIFvPatchScalarField
(
    const energyJumpAMIFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedJumpAMIFvPatchField<scalar>(ptf, p, iF, mapper)
{}


Foam::energyJumpAMIFvPatchScalarField::energyJumpAMIFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedJumpAMIFvPatchField<scalar>(p, iF)
{
    FatalErrorInFunction
        << "Invalid assignment of the 'energyJumpAMI' patch field type "
        << "directly as a field boundary condition. Please assign the "
        << "'fixedJumpAMI' patch field type instead."
        << exit(FatalError);
}


Foam::energyJumpAMIFvPatchScalarField::energyJumpAMIFvPatchScalarField
(
    const energyJumpAMIFvPatchScalarField& ptf
)
:
    fixedJumpAMIFvPatchField<scalar>(ptf)
{}


Foam::energyJumpAMIFvPatchScalarField::energyJumpAMIFvPatchScalarField
(
    const energyJumpAMIFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpAMIFvPatchField<scalar>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::energyJumpAMIFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (this->cyclicAMIPatch().owner())
    {
        const basicThermo& thermo = basicThermo::lookupThermo(*this);
        const label patchi = patch().index();
        const fixedJumpAMIFvPatchScalarField& TbPatch =
            refCast<const fixedJumpAMIFvPatchScalarField>
            (
                thermo.T().boundaryField()[patchi]
            );

        fixedJumpAMIFvPatchScalarField& Tbp =
            const_cast<fixedJumpAMIFvPatchScalarField&>(TbPatch);

        // force update of jump
        Tbp.updateCoeffs();

        const labelUList& faceCells = this->patch().faceCells();
        const scalarField& pp = thermo.p().boundaryField()[patchi];
        jump_ = thermo.he(pp, Tbp.jump(), faceCells);
    }

    fixedJumpAMIFvPatchField<scalar>::updateCoeffs();
}


void Foam::energyJumpAMIFvPatchScalarField::write(Ostream& os) const
{
    fixedJumpAMIFvPatchField<scalar>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchScalarField,
       energyJumpAMIFvPatchScalarField
   );
}

// ************************************************************************* //
