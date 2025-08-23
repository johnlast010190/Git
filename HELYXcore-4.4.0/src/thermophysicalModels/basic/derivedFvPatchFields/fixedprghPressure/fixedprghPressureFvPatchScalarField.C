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
    (c) 2020-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "derivedFvPatchFields/fixedprghPressure/fixedprghPressureFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "basicThermo/basicThermo.H"


// * * * * * * * * * * * Protected member functions  * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::fixedprghPressureFvPatchScalarField::patchInternalRho()
{
    basicThermo& thermo =
        db().lookupObjectRef<basicThermo>(basicThermo::matDictName);

    const fvPatch& fvP = patch();

    tmp<scalarField> prho
    (
        new scalarField
        (
            thermo.distinctBuoyancy()
          ? thermo.buoyantRho()().boundaryField()
            [
                fvP.index()
            ].patchInternalField()
          : thermo.rho()().boundaryField()[fvP.index()].patchInternalField()
        )
    );

    return prho;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedprghPressureFvPatchScalarField::fixedprghPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    p_rgh_(p.size(), 0.0),
    hRef_(0.0),
    hRefSpecified_(false)
{}


Foam::fixedprghPressureFvPatchScalarField::fixedprghPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    p_rgh_("p_rgh", dict, p.size()),
    hRef_(dict.lookupOrDefault<scalar>("hRef", 0.0)),
    hRefSpecified_(dict.found("hRef"))
{}


Foam::fixedprghPressureFvPatchScalarField::fixedprghPressureFvPatchScalarField
(
    const fixedprghPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    p_rgh_(mapper(ptf.p_rgh_)),
    hRef_(ptf.hRef_),
    hRefSpecified_(ptf.hRefSpecified_)
{}


Foam::fixedprghPressureFvPatchScalarField::fixedprghPressureFvPatchScalarField
(
    const fixedprghPressureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    p_rgh_(tppsf.p_rgh_),
    hRef_(tppsf.hRef_),
    hRefSpecified_(tppsf.hRefSpecified_)
{}


Foam::fixedprghPressureFvPatchScalarField::fixedprghPressureFvPatchScalarField
(
    const fixedprghPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    p_rgh_(tppsf.p_rgh_),
    hRef_(tppsf.hRef_),
    hRefSpecified_(tppsf.hRefSpecified_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedprghPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    m(p_rgh_, p_rgh_);
}


void Foam::fixedprghPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const fixedprghPressureFvPatchScalarField& tiptf =
        refCast<const fixedprghPressureFvPatchScalarField>(ptf);

    p_rgh_.rmap(tiptf.p_rgh_, addr);
}


void Foam::fixedprghPressureFvPatchScalarField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchScalarField::autoMapGIB(mapper);
    mapper.map(p_rgh_, gAverage(p_rgh_));
}


void Foam::fixedprghPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalarField pp(p_rgh_);

    // Do not add if solving for potential flow, as gravity
    // is not accounted for
    if (this->internalField().name() == "pPot")
    {
        return;
    }

    dimensionedVector g("g", dimAcceleration, vector::zero);
    const uniformDimensionedVectorField* gPtr =
        this->db().lookupObjectPtr<uniformDimensionedVectorField>("g");

    const fvPatch& fvP = this->patch();

    if (gPtr)
    {
        g = *gPtr;
    }

    // Points adjacent to the cell centres, rather than face centres
    pointField patchFacePoints
    (
        fvP.patchInternalField(this->internalField().mesh().C())
      + fvP.delta()
    );

    scalarField gh
    (
        g.value() & patchFacePoints
    );

    // Use specified hRef or solver hRef (if present)
    dimensionedScalar hRef("hRef", dimLength, 0.0);
    if (hRefSpecified_)
    {
        hRef = dimensionedScalar("hRef", dimLength, hRef_);
    }
    else
    {
        const uniformDimensionedScalarField* hRefPtr =
            db().lookupObjectPtr<uniformDimensionedScalarField>("hRef");
        if (hRefPtr)
        {
            hRef = *hRefPtr;
        }
    }

    dimensionedScalar ghRef
    (
        mag(g.value()) > SMALL
        ? g & (cmptMag(g.value())/mag(g.value()))*hRef
        : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    );

    scalarField rhop(patchInternalRho());
    pp += rhop*(gh - ghRef.value());

    forceAssign(pp);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::fixedprghPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    referenceFrameFvPatch::write(os);
    p_rgh_.writeEntry("p_rgh", os);
    if (hRefSpecified_)
    {
        os.writeEntry("hRef", hRef_);
    }
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedprghPressureFvPatchScalarField
    );
}

// ************************************************************************* //
