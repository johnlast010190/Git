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
    (c) 2022-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "derivedFvPatchFields/coupledTotalPressure/coupledTotalPressureFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "addStaticHead/staticHead.H"
#include "cfdTools/general/solutionControl/helyxCoupledControl/helyxCoupledControl.H"

// * * * * * * * * * * * * * Private Functions   * * * * * * * * * * * * * * //

void Foam::coupledTotalPressureFvPatchScalarField::computeOtherSources()
{
    if (staticHead_.active())
    {
        pStaticHead_ = staticHead_.computeStaticHead();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledTotalPressureFvPatchScalarField::
coupledTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    basePressureFvPatchScalarField(p, iF),
    pStaticHead_(),
    staticHead_(*this),
    freeStreamVelocity_(nullptr)
{}


Foam::coupledTotalPressureFvPatchScalarField::
coupledTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    basePressureFvPatchScalarField(p, iF, dict),
    pStaticHead_(),
    staticHead_(*this, p, dict),
    freeStreamVelocity_
    (
        dict.found("freeStreamVelocity")
        ? new point(dict.lookup("freeStreamVelocity"))
        : nullptr
    )
{
    if (staticHead_.active())
    {
        if (dict.found("pStaticHead"))
        {
            pStaticHead_ = scalarField("pStaticHead", dict, p.size());
        }
        else
        {
            pStaticHead_ = scalarField(p.size(), 0);
        }
    }
}


Foam::coupledTotalPressureFvPatchScalarField::
coupledTotalPressureFvPatchScalarField
(
    const coupledTotalPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    basePressureFvPatchScalarField(ptf, p, iF, mapper),
    pStaticHead_(mapper(ptf.pStaticHead_)),
    staticHead_(*this, ptf.staticHead_, mapper),

    freeStreamVelocity_
    (
         ptf.freeStreamVelocity_.valid()
       ? new vector(ptf.freeStreamVelocity_())
       : nullptr
    )
{}


Foam::coupledTotalPressureFvPatchScalarField::
coupledTotalPressureFvPatchScalarField
(
    const coupledTotalPressureFvPatchScalarField& tppsf
)
:
    basePressureFvPatchScalarField(tppsf),
    pStaticHead_(tppsf.pStaticHead_),
    staticHead_(*this, tppsf.staticHead_),
    freeStreamVelocity_
    (
         tppsf.freeStreamVelocity_.valid()
       ? new vector(tppsf.freeStreamVelocity_())
       : nullptr
    )
{}


Foam::coupledTotalPressureFvPatchScalarField::
coupledTotalPressureFvPatchScalarField
(
    const coupledTotalPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    basePressureFvPatchScalarField(tppsf, iF),
    pStaticHead_(tppsf.pStaticHead_),
    staticHead_(*this, tppsf.staticHead_),
    freeStreamVelocity_
    (
         tppsf.freeStreamVelocity_.valid()
       ? new vector(tppsf.freeStreamVelocity_())
       : nullptr
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Mapping functions

//- Map (and resize as needed) from self given a mapping object
void Foam::coupledTotalPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    basePressureFvPatchScalarField::autoMap(m);

    if (staticHead_.active())
    {
        m(pStaticHead_, pStaticHead_);
        staticHead_.autoMap(m);
    }
}


//- Reverse map the given fvPatchField onto this fvPatchField
void Foam::coupledTotalPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    basePressureFvPatchScalarField::rmap(ptf, addr);

    const coupledTotalPressureFvPatchScalarField& tiptf =
        refCast<const coupledTotalPressureFvPatchScalarField>(ptf);

    if (staticHead_.active())
    {
        pStaticHead_.rmap(tiptf.pStaticHead_, addr);
        staticHead_.rmap(tiptf.staticHead_, addr);
    }
}


//- Map (and resize as needed) from self given a mapping object
void Foam::coupledTotalPressureFvPatchScalarField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    basePressureFvPatchScalarField::autoMapGIB(mapper);
    if (staticHead_.active())
    {
        mapper.map(pStaticHead_, scalar(0));
        staticHead_.autoMapGIB(mapper);
    }
}



Foam::tmp<Foam::scalarField>
Foam::coupledTotalPressureFvPatchScalarField::C0Field() const
{
    tmp<scalarField> tc0(new scalarField(this->p0()));
    scalarField& c0 = tc0.ref();

    if (staticHead_.active())
    {
        c0 += pStaticHead_;
    }

    if (freeStreamVelocity_.valid())
    {
        scalarField c0fs(c0.size(), Zero);

        const vectorField UpIn
        (
            this->db().lookupObject<volVectorField>(UName()).
                boundaryField()[this->patch().index()]
        );
        scalar magSqUfreeSteam = magSqr(freeStreamVelocity_());

        // Urel^2 = (U_i - Ufree_i)^2 = Ui^2-2*U_i&Ufree_i + Ufree_i^2;
        // Ui^2 is added implicitly in c2 and c3
        c0fs = neg(phiP())*((freeStreamVelocity_()&UpIn) - 0.5*magSqUfreeSteam);
        if (internalField().dimensions() == dimPressure)
        {
            c0fs *= this->lookupPatchField<volScalarField, scalar>(rhoName());
        }

        c0 += c0fs;
    }
    return tc0;
}


Foam::tmp<Foam::scalarField>
Foam::coupledTotalPressureFvPatchScalarField::C1Field() const
{
    return tmp<scalarField>
    (
        new scalarField(this->size(), Zero)
    );
}


Foam::tmp<Foam::scalarField>
Foam::coupledTotalPressureFvPatchScalarField::C2Field() const
{
    tmp<scalarField> tc2
    (
        new scalarField(this->size(), -0.5)
    );
    scalarField& c2 = tc2.ref();

    scalarField relPhip(phiP());
    if (!this->db().foundObject<helyxCoupledControl>(solutionControl::typeName))
    {
        if (internalField().dimensions() == dimPressure)
        {
            const fvPatchField<scalar>& rho =
                patch().lookupPatchFieldInDb<volScalarField, scalar>
                (
                    this->db(), staticHead_.getRhoName()
                );
            scalarField phiRef(this->size(), Zero);
            makeRelative(phiRef);
            phiRef *= rho;
            relPhip += phiRef;
        }
        else
        {
            makeRelative(relPhip);
        }
    }
    c2 *= neg(relPhip);

    if (internalField().dimensions() == dimPressure)
    {
        const fvPatchField<scalar>& rho =
            patch().lookupPatchFieldInDb<volScalarField, scalar>
            (
                this->db(), staticHead_.getRhoName()
            );
        c2 *= rho;
    }

    return tc2;
}


Foam::tmp<Foam::scalarField>
Foam::coupledTotalPressureFvPatchScalarField::C3Field() const
{
    tmp<scalarField> tc3
    (
        new scalarField(this->size(), -0.5)
    );
    scalarField& c3 = tc3.ref();
    scalarField relPhip(phiP());

    if (!this->db().foundObject<helyxCoupledControl>(solutionControl::typeName))
    {
        if (internalField().dimensions() == dimPressure)
        {
            const fvPatchField<scalar>& rho =
                patch().lookupPatchFieldInDb<volScalarField, scalar>
                (
                    this->db(), staticHead_.getRhoName()
                );
            scalarField phiRef(this->size(), Zero);
            makeRelative(phiRef);
            phiRef *= rho;
            relPhip += phiRef;
        }
        else
        {
            makeRelative(relPhip);
        }
    }
    c3 *= neg(relPhip);

    if (internalField().dimensions() == dimPressure)
    {
        const fvPatchField<scalar>& rho =
            patch().lookupPatchFieldInDb<volScalarField, scalar>
            (
                this->db(), staticHead_.getRhoName()
            );
        c3 *= rho;
    }

    return tc3;
}


void Foam::coupledTotalPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    referenceFrameFvPatch::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName());
    p0().writeEntry("p0", os);
    staticHead_.write(os);
    if (staticHead_.active())
    {
        pStaticHead_.writeEntry("pStaticHead", os);
    }
    if (freeStreamVelocity_.valid())
    {
        os.writeEntry("freeStreamVelocity", freeStreamVelocity_());
    }
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        coupledTotalPressureFvPatchScalarField
    );
}

// ************************************************************************* //
