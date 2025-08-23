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

#include "turbulentFluidThermoModels/derivedFvPatchFields/wallFunctions/alphatWallFunctions/alphatPhaseChangeWallFunction/alphatPhaseChangeWallFunctionFvPatchScalarField.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(alphatPhaseChangeWallFunctionFvPatchScalarField,0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatPhaseChangeWallFunctionFvPatchScalarField::
alphatPhaseChangeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    dmdt_(p.size(), 0),
    mDotL_(p.size(), 0)
{}


alphatPhaseChangeWallFunctionFvPatchScalarField::
alphatPhaseChangeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    dmdt_(p.size(), 0),
    mDotL_(p.size(), 0)
{
    if (dict.found("dmdt"))
    {
        dmdt_ = scalarField("dmdt", dict, p.size());
    }

    if (dict.found("mDotL"))
    {
        dmdt_ = scalarField("mDotL", dict, p.size());
    }
}


alphatPhaseChangeWallFunctionFvPatchScalarField::
alphatPhaseChangeWallFunctionFvPatchScalarField
(
    const alphatPhaseChangeWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    dmdt_(mapper(ptf.dmdt_)),
    mDotL_(mapper(ptf.mDotL_))
{}


alphatPhaseChangeWallFunctionFvPatchScalarField::
alphatPhaseChangeWallFunctionFvPatchScalarField
(
    const alphatPhaseChangeWallFunctionFvPatchScalarField& awfpsf
)
:
    fixedValueFvPatchScalarField(awfpsf),
    dmdt_(awfpsf.dmdt_),
    mDotL_(awfpsf.mDotL_)
{}


alphatPhaseChangeWallFunctionFvPatchScalarField::
alphatPhaseChangeWallFunctionFvPatchScalarField
(
    const alphatPhaseChangeWallFunctionFvPatchScalarField& awfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(awfpsf, iF),
    dmdt_(awfpsf.dmdt_),
    mDotL_(awfpsf.mDotL_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void alphatPhaseChangeWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);

    m(dmdt_, dmdt_);
    m(mDotL_, mDotL_);
}


void alphatPhaseChangeWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const alphatPhaseChangeWallFunctionFvPatchScalarField& tiptf =
        refCast<const alphatPhaseChangeWallFunctionFvPatchScalarField>(ptf);

    dmdt_.rmap(tiptf.dmdt_, addr);
    mDotL_.rmap(tiptf.mDotL_, addr);
}


void alphatPhaseChangeWallFunctionFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    fixedValueFvPatchScalarField::reset(ptf);

    const alphatPhaseChangeWallFunctionFvPatchScalarField& tiptf =
        refCast<const alphatPhaseChangeWallFunctionFvPatchScalarField>(ptf);

    dmdt_.reset(tiptf.dmdt_);
    mDotL_.reset(tiptf.mDotL_);
}


void alphatPhaseChangeWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    dmdt_.writeEntry("dmdt", os);
    mDotL_.writeEntry("mDotL", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
