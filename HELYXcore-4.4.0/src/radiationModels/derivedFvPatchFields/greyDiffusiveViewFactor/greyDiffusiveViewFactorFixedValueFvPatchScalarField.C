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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "derivedFvPatchFields/greyDiffusiveViewFactor/greyDiffusiveViewFactorFixedValueFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "radiationModels/radiationModel/radiationModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::greyDiffusiveViewFactorFixedValueFvPatchScalarField::
greyDiffusiveViewFactorFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    qro_(),
    solarLoad_(false)
{}


Foam::radiationModels::greyDiffusiveViewFactorFixedValueFvPatchScalarField::
greyDiffusiveViewFactorFixedValueFvPatchScalarField
(
    const greyDiffusiveViewFactorFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    qro_(mapper(ptf.qro_)),
    solarLoad_(ptf.solarLoad_)
{}


Foam::radiationModels::greyDiffusiveViewFactorFixedValueFvPatchScalarField::
greyDiffusiveViewFactorFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    qro_("qro", dict, p.size()),
    solarLoad_(dict.lookupOrDefault<bool>("solarLoad", false))
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
        fvPatchScalarField::operator=(0.0);
    }
}


Foam::radiationModels::greyDiffusiveViewFactorFixedValueFvPatchScalarField::
greyDiffusiveViewFactorFixedValueFvPatchScalarField
(
    const greyDiffusiveViewFactorFixedValueFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    qro_(ptf.qro_),
    solarLoad_(ptf.solarLoad_)
{}


Foam::radiationModels::greyDiffusiveViewFactorFixedValueFvPatchScalarField::
greyDiffusiveViewFactorFixedValueFvPatchScalarField
(
    const greyDiffusiveViewFactorFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    qro_(ptf.qro_),
    solarLoad_(ptf.solarLoad_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiationModels::greyDiffusiveViewFactorFixedValueFvPatchScalarField::
autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    m(qro_, qro_);
}


void Foam::radiationModels::greyDiffusiveViewFactorFixedValueFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const greyDiffusiveViewFactorFixedValueFvPatchScalarField& mrptf =
        refCast<const greyDiffusiveViewFactorFixedValueFvPatchScalarField>(ptf);

    qro_.rmap(mrptf.qro_, addr);
}


void Foam::radiationModels::greyDiffusiveViewFactorFixedValueFvPatchScalarField::
reset
(
    const fvPatchScalarField& ptf
)
{
    fixedValueFvPatchScalarField::reset(ptf);

    const greyDiffusiveViewFactorFixedValueFvPatchScalarField& mrptf =
        refCast<const greyDiffusiveViewFactorFixedValueFvPatchScalarField>(ptf);

    qro_.reset(mrptf.qro_);
}


void Foam::radiationModels::greyDiffusiveViewFactorFixedValueFvPatchScalarField::
autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchScalarField::autoMapGIB(mapper);
    mapper.map(qro_, scalar(0));
}


void Foam::radiationModels::greyDiffusiveViewFactorFixedValueFvPatchScalarField::
updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (debug)
    {
        scalar Q = gSum((*this)*patch().magSf());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << " heat transfer rate:" << Q
            << " wall radiative heat flux "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}


Foam::tmp<Foam::scalarField> Foam::radiationModels::
greyDiffusiveViewFactorFixedValueFvPatchScalarField::qro() const
{
    tmp<scalarField> tqrt(new scalarField(qro_));

    if (solarLoad_)
    {
        const radiationModel& radiation =
            db().lookupObject<radiationModel>("radiationProperties");

        tqrt.ref() += patch().lookupPatchField<volScalarField, scalar>
        (
            radiation.externalRadHeatFieldName_
        );
    }

    return tqrt;
}


void Foam::radiationModels::greyDiffusiveViewFactorFixedValueFvPatchScalarField::
write
(
    Ostream& os
) const
{
    fixedValueFvPatchScalarField::write(os);
    qro_.writeEntry("qro", os);
    os.writeEntry("solarLoad", solarLoad_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
    makePatchTypeField
    (
        fvPatchScalarField,
        greyDiffusiveViewFactorFixedValueFvPatchScalarField
    );
}
}


// ************************************************************************* //
