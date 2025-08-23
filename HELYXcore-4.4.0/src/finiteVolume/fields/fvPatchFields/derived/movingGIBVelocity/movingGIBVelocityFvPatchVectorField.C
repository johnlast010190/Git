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
    (c) 2018-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/movingGIBVelocity/movingGIBVelocityFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "finiteVolume/fvc/fvcMeshPhi.H"
#include "fvMesh/fvMeshGIBChangers/fvMeshGIBChanger/fvMeshGIBChanger.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingGIBVelocityFvPatchVectorField::
movingGIBVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    UName_("U")
{}


Foam::movingGIBVelocityFvPatchVectorField::
movingGIBVelocityFvPatchVectorField
(
    const movingGIBVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    UName_(ptf.UName_)
{}


Foam::movingGIBVelocityFvPatchVectorField::
movingGIBVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U"))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::movingGIBVelocityFvPatchVectorField::
movingGIBVelocityFvPatchVectorField
(
    const movingGIBVelocityFvPatchVectorField& mwvpvf
)
:
    fixedValueFvPatchVectorField(mwvpvf),
    UName_(mwvpvf.UName_)
{}


Foam::movingGIBVelocityFvPatchVectorField::
movingGIBVelocityFvPatchVectorField
(
    const movingGIBVelocityFvPatchVectorField& mwvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(mwvpvf, iF),
    UName_(mwvpvf.UName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::movingGIBVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = internalField().mesh();

    if (mesh.changing())
    {
        if (mesh.hasChangers())
        {
            vectorField::operator=
            (
                mesh.GIBChanger().velocityCorrect(this->patch().Cf())
            );
        }
        else
        {
            vectorField::operator=(mesh.velocityCorrect(this->patch().Cf()));
        }
    }
    else
    {
        vectorField::operator=(vector::zero);
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::movingGIBVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        movingGIBVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
