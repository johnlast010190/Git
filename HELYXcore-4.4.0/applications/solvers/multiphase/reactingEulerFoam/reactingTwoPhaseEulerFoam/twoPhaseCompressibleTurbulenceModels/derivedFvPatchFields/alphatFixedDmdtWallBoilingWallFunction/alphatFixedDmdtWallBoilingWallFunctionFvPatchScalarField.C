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
    (c) 2015-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::
alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(p, iF),
    relax_(1.0),
    fixedDmdt_(0.0),
    L_(0.0)
{
    checkType();
}


alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::
alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(p, iF, dict),
    relax_(dict.lookupOrDefault<scalar>("relax", 1.0)),
    fixedDmdt_(dict.lookupOrDefault<scalar>("fixedDmdt", 0.0)),
    L_(dict.lookupOrDefault<scalar>("L", 0.0))
{}


alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::
alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
(
    const alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField
    (
        psf,
        p,
        iF,
        mapper
    ),
    fixedDmdt_(psf.fixedDmdt_),
    L_(psf.L_)
{}


alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::
alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
(
    const alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField& psf
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(psf),
    relax_(psf.relax_),
    fixedDmdt_(psf.fixedDmdt_),
    L_(psf.L_)
{}

alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::
alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
(
    const alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(psf, iF),
    relax_(psf.relax_),
    fixedDmdt_(psf.fixedDmdt_),
    L_(psf.L_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    dmdt_ = (1 - relax_)*dmdt_ + relax_*fixedDmdt_;
    mDotL_ = dmdt_*L_;

    forceAssign(calcAlphat(*this));

    fixedValueFvPatchScalarField::updateCoeffs();
}


void alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry("relax", relax_);
    os.writeEntry("fixedDmdt", fixedDmdt_);
    os.writeEntry("L", L_);
    dmdt_.writeEntry("dmdt", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
);

// Backward compatibility
addSpecialNamedToPatchFieldRunTimeSelection
(
    fvPatchScalarField,
    alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField,
    compressible::alphatFixedDmdtWallBoilingWallFunction,
    compressible__alphatFixedDmdtWallBoilingWallFunction
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
