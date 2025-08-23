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
    (c) 2023-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/movingWallSlipVelocity/movingWallSlipVelocityFvPatchVectorField.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "finiteVolume/fvc/fvcMeshPhi.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingWallSlipVelocityFvPatchVectorField::
movingWallSlipVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF)
{
    const vectorField n(p.nf());

    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = sqr(n);
}


Foam::movingWallSlipVelocityFvPatchVectorField::
movingWallSlipVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

    const vectorField n(p.nf());

    refValue() = sqr(n) & *this;
    refGrad() = Zero;
    valueFraction() = sqr(n);
}


Foam::movingWallSlipVelocityFvPatchVectorField::
movingWallSlipVelocityFvPatchVectorField
(
    const movingWallSlipVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::movingWallSlipVelocityFvPatchVectorField::
movingWallSlipVelocityFvPatchVectorField
(
    const movingWallSlipVelocityFvPatchVectorField& mwsvpvf
)
:
    directionMixedFvPatchVectorField(mwsvpvf)
{}


Foam::movingWallSlipVelocityFvPatchVectorField::
movingWallSlipVelocityFvPatchVectorField
(
    const movingWallSlipVelocityFvPatchVectorField& mwsvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(mwsvpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::movingWallSlipVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = internalField().mesh();

    if (mesh.moving())
    {
        const fvPatch& p = patch();

        const volVectorField& U =
            static_cast<const volVectorField&>(internalField());

        const scalarField phip(fvc::meshPhi(U, p.index()));

        const vectorField n(p.nf());
        const scalarField& magSf = p.magSf();

        tmp<scalarField> Un = phip/(magSf + VSMALL);

        refValue() = n*Un;
        refGrad() = Zero;
        valueFraction() = sqr(n);
    }

    directionMixedFvPatchVectorField::updateCoeffs();
    directionMixedFvPatchVectorField::evaluate();
}


void Foam::movingWallSlipVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        movingWallSlipVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
