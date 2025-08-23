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
#include "derivedFvPatchFields/energyJump/energyJump/energyJumpFvPatchScalarField.H"
#include "fields/fvPatchFields/derived/fixedJump/fixedJumpFvPatchFields.H"
#include "basicThermo/basicThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::energyJumpFvPatchScalarField::energyJumpFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpFvPatchField<scalar>(p, iF)
{}


Foam::energyJumpFvPatchScalarField::energyJumpFvPatchScalarField
(
    const energyJumpFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedJumpFvPatchField<scalar>(ptf, p, iF, mapper)
{}


Foam::energyJumpFvPatchScalarField::energyJumpFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedJumpFvPatchField<scalar>(p, iF)
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        evaluate(Pstream::commsTypes::blocking);
    }
}


Foam::energyJumpFvPatchScalarField::energyJumpFvPatchScalarField
(
    const energyJumpFvPatchScalarField& ptf
)
:
    fixedJumpFvPatchField<scalar>(ptf)
{}


Foam::energyJumpFvPatchScalarField::energyJumpFvPatchScalarField
(
    const energyJumpFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpFvPatchField<scalar>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::energyJumpFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (this->cyclicPatch().owner())
    {
        const basicThermo& thermo = basicThermo::lookupThermo(*this);
        const label patchi = patch().index();
        const scalarField& pp = thermo.p().boundaryField()[patchi];

        const fixedJumpFvPatchScalarField& TbPatch =
            refCast<const fixedJumpFvPatchScalarField>
            (
                thermo.T().boundaryField()[patchi]
            );

        fixedJumpFvPatchScalarField& Tbp =
            const_cast<fixedJumpFvPatchScalarField&>(TbPatch);

        // force update of jump
        Tbp.updateCoeffs();

        jump_ = thermo.he(pp, Tbp.jump(), this->patch().faceCells());
    }

    fixedJumpFvPatchField<scalar>::updateCoeffs();
}


void Foam::energyJumpFvPatchScalarField::write(Ostream& os) const
{
    fixedJumpFvPatchField<scalar>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchScalarField,
       energyJumpFvPatchScalarField
   );
}

// ************************************************************************* //
