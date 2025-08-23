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
    (c) 2018-2023 OpenFOAM Foundation
    (c) 2023-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/freestreamVelocity/freestreamVelocityFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::freestreamVelocityFvPatchVectorField::freestreamVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF)
{}


Foam::freestreamVelocityFvPatchVectorField::freestreamVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF, dict, false)
{
    referenceFrameFvPatch<vector>::read(dict);
    inletFlux_ = true;
    this->patchType() = dict.lookupOrDefault<word>("patchType", word::null);

    if (dict.found("freestreamValueF1"))
    {
        freestreamValue_.reset(Function1<vector>::New("freestreamValueF1", dict));
        freestreamValue() =
            freestreamValue_->value(this->db().time().timeOutputValue());
    }
    else
    {
        freestreamValue() = vectorField("freestreamValue", dict, p.size());
    }

    if (dict.found("value"))
    {
        fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        fvPatchVectorField::operator=(freestreamValue());
    }

    refGrad() = Zero;
    if (dict.found("valueFraction"))
    {
        valueFraction()  = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        valueFraction() = 1;
    }
}


Foam::freestreamVelocityFvPatchVectorField::freestreamVelocityFvPatchVectorField
(
    const freestreamVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    freestreamValue_(ptf.freestreamValue_, false)
{}


Foam::freestreamVelocityFvPatchVectorField::freestreamVelocityFvPatchVectorField
(
    const freestreamVelocityFvPatchVectorField& ptf
)
:
    mixedFvPatchVectorField(ptf),
    freestreamValue_(ptf.freestreamValue_, false)
{}


Foam::freestreamVelocityFvPatchVectorField::freestreamVelocityFvPatchVectorField
(
    const freestreamVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(ptf, iF),
    freestreamValue_(ptf.freestreamValue_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::freestreamVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchVectorField::autoMap(m);

    // Override
    if (freestreamValue_.valid())
    {
        freestreamValue() =
            freestreamValue_->value(this->db().time().timeOutputValue());
    }
}


void Foam::freestreamVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    mixedFvPatchVectorField::rmap(ptf, addr);

    // Override
    if (freestreamValue_.valid())
    {
        freestreamValue() =
            freestreamValue_->value(this->db().time().timeOutputValue());
    }
}


void Foam::freestreamVelocityFvPatchVectorField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    mixedFvPatchVectorField::autoMapGIB(mapper);

    // Override
    if (freestreamValue_.valid())
    {
        freestreamValue() =
            freestreamValue_->value(this->db().time().timeOutputValue());
    }
}


void Foam::freestreamVelocityFvPatchVectorField::reset
(
    const fvPatchVectorField& ptf
)
{
    mixedFvPatchVectorField::reset(ptf);

    // Override
    if (freestreamValue_.valid())
    {
        freestreamValue() =
            freestreamValue_->value(this->db().time().timeOutputValue());
    }
}


void Foam::freestreamVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Update the uniform value as a function of time
    if (freestreamValue_.valid())
    {
        const scalar t = this->db().time().timeOutputValue();
        freestreamValue() = freestreamValue_->value(t);
    }

    vectorField Up(0.5*(patchInternalField() + *this));

    //- here refFrames can be a bit complicated
    //  Three possible configurations
    //  1. referenceFrame not specified - classic approach
    //  2. referenceFrame specified and the whole domain in MRF. Here, we want
    //     the angle to be computed based on Urel and not U, so we add the frame
    //     velocity to the calculation of angle.
    //  3. referenceFrame specified but used to calculate the freestream
    //     velocity on a frame but without MRF. Here U = Urel. So we need the
    //     frameVelocity to be accounted in the U evaluation.

    if (coorFramePtr())
    {
        coorFramePtr()->attachPatch(this->patch().index());
        if (coorFramePtr()->isAttachToMRF(this->patch().index()))
        {
            vectorField Uframe(this->patch().size(), Zero);
            coorFramePtr()->inletFluxVelocity
            (
                Uframe,
                this->patch().index(),
                this->patch().boundaryMesh().mesh(),
                this->internalField().name()
            );
            Up -= Uframe;
        }
    }

    const scalarField magUp(mag(Up));
    const vectorField nf(patch().nf());

    scalarField& vf = valueFraction();

    forAll(vf, i)
    {
        if (magUp[i] > VSMALL)
        {
            vf[i] = 0.5 - 0.5*(Up[i] & nf[i])/magUp[i];
        }
        else
        {
            vf[i] = 0.5;
        }
    }

    mixedFvPatchField<vector>::updateCoeffs();
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::freestreamVelocityFvPatchVectorField::getFrameVelocity() const
{
    tmp<vectorField> tframe(new vectorField(patch().size(), Zero));
    if (coorFramePtr())
    {
        if (!coorFramePtr()->isAttachToMRF(this->patch().index()))
        {
            tframe = mixedFvPatchField<vector>::getFrameVelocity();
        }
    }
    return tframe;
}


void Foam::freestreamVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    freestreamValue().writeEntry("freestreamValue", os);
    if (freestreamValue_.valid())
    {
        freestreamValue_->writeData(os);
    }
    valueFraction().writeEntry("valueFraction", os);
    referenceFrameFvPatch<vector>::write(os);
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        freestreamVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
