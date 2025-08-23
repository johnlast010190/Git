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
    (c) 2019-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/wallVelocity/wallVelocityFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "finiteVolume/fvc/fvcMeshPhi.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallVelocityFvPatchVectorField::wallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::wallVelocityFvPatchVectorField::wallVelocityFvPatchVectorField
(
    const wallVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::wallVelocityFvPatchVectorField::wallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false)
{
    if (dict.found("value"))
    {
        forceAssign(vectorField("value", dict, p.size()));
    }
    else
    {
        wallVelocityFvPatchVectorField::updateCoeffs();
    }
}


Foam::wallVelocityFvPatchVectorField::wallVelocityFvPatchVectorField
(
    const wallVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf)
{}


Foam::wallVelocityFvPatchVectorField::wallVelocityFvPatchVectorField
(
    const wallVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wallVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    vectorField Up(inputValue());
    frameFieldUpdate(Up, false);

    fixedValueFvPatchVectorField::updateCoeffs();
}


Foam::tmp<Foam::CoeffField<Foam::vector>>
Foam::wallVelocityFvPatchVectorField::gradientInternalBCoeffs() const
{
    return fixedValueFvPatchField<vector>::gradientInternalBCoeffs();
}


Foam::tmp<Foam::vectorField>
Foam::wallVelocityFvPatchVectorField::gradientBoundaryCoeffs() const
{
    return fixedValueFvPatchField<vector>::gradientBoundaryCoeffs();
}


Foam::tmp<Foam::vectorField>
Foam::wallVelocityFvPatchVectorField::gradientBoundaryBCoeffs() const
{
    return this->gradientBoundaryCoeffs();
}


Foam::tmp<Foam::vectorField>
Foam::wallVelocityFvPatchVectorField::snGrad() const
{
    tmp<vectorField> nHat = this->patch().nf();
    vectorField dfield(*this - this->patchInternalField());
    dfield -= (dfield & nHat())*nHat();
    return (dfield*this->patch().deltaCoeffs());
}


void Foam::wallVelocityFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField(fvPatchVectorField, wallVelocityFvPatchVectorField);
}

// ************************************************************************* //
