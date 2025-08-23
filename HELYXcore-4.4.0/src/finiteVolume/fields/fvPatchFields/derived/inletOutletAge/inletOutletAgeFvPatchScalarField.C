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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/inletOutletAge/inletOutletAgeFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inletOutletAgeFvPatchScalarField::inletOutletAgeFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


Foam::inletOutletAgeFvPatchScalarField::inletOutletAgeFvPatchScalarField
(
    const inletOutletAgeFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::inletOutletAgeFvPatchScalarField::inletOutletAgeFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF, dict)
{
    this->refValue() = this->patch().boundaryMesh().mesh().time().value();
    fvPatchField<scalar>::forceAssign(this->refValue());
}


Foam::inletOutletAgeFvPatchScalarField::inletOutletAgeFvPatchScalarField
(
    const inletOutletAgeFvPatchScalarField& ptf
)
:
    inletOutletFvPatchScalarField(ptf)
{}


Foam::inletOutletAgeFvPatchScalarField::inletOutletAgeFvPatchScalarField
(
    const inletOutletAgeFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::inletOutletAgeFvPatchScalarField::updateCoeffs()
{
    this->refValue() = this->patch().boundaryMesh().mesh().time().value();

    inletOutletFvPatchScalarField::updateCoeffs();
}


void Foam::inletOutletAgeFvPatchScalarField::write(Ostream& os) const
{
    inletOutletFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::inletOutletAgeFvPatchScalarField::operator=
(
    const fvPatchField<scalar>& ptf
)
{
    fvPatchField<scalar>::forceAssign
    (
        this->valueFraction()*this->patch().boundaryMesh().mesh().time().value()
        + (1 - this->valueFraction())*ptf
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        inletOutletAgeFvPatchScalarField
    );
}

// ************************************************************************* //
