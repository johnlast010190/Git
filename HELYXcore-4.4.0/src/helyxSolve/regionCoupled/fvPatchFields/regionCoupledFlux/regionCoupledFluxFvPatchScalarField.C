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
    (c) 2012-2019 OpenFOAM Foundation
    (c) 2010-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "regionCoupledFluxFvPatchScalarField.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"
#include "db/Time/Time.H"
#include "fields/volFields/volFields.H"
#include "fvMatrices/fvMatrices.H"


// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

regionCoupledFluxFvPatchScalarField::
regionCoupledFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCoupledFvPatchField<scalar>(p, iF),
    faceBCoeff_(p.size(), 0),
    faceICoeff_(p.size(), 1),
    faceCorr_(p.size(), 0),
    contactResistance_(0),
    contactInsulance_(0),
    Qsum_(p.size(), 0.0),
    QsumNbr_(p.size(), 0.0),
    qURF_(0.3)
{}


regionCoupledFluxFvPatchScalarField::
regionCoupledFluxFvPatchScalarField
(
    const regionCoupledFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    regionCoupledFvPatchField<scalar>(ptf, p, iF, mapper),
    faceBCoeff_(mapper(ptf.faceBCoeff_)),
    faceICoeff_(mapper(ptf.faceICoeff_)),
    faceCorr_(mapper(ptf.faceCorr_)),
    contactResistance_(ptf.contactResistance_),
    conductivityLayers_(ptf.conductivityLayers_),
    thicknessLayers_(ptf.thicknessLayers_),
    contactInsulance_(ptf.contactInsulance_),
    Qsum_(mapper(ptf.Qsum_)),
    QsumNbr_(mapper(ptf.QsumNbr_)),
    qURF_(ptf.qURF_)
{}


regionCoupledFluxFvPatchScalarField::
regionCoupledFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    regionCoupledFvPatchField<scalar>(p, iF, dict),
    faceBCoeff_
    (
        dict.found("faceBCoeff")
      ? std::move(scalarField(word("faceBCoeff"), dict, p.size()))
      : std::move(scalarField(p.size(), 0.0))
    ),
    faceICoeff_
    (
        dict.found("faceICoeff")
      ? std::move(scalarField(word("faceICoeff"), dict, p.size()))
      : std::move(scalarField(p.size(), 1.0))
    ),
    faceCorr_
    (
        dict.found("faceCorr")
      ? std::move(scalarField(word("faceCorr"), dict, p.size()))
      : std::move(scalarField(p.size(), 0.0))
    ),
    contactResistance_(scalar(0)),
    contactInsulance_(scalar(0)),
    Qsum_
    (
        dict.found("Qsum")
      ? std::move(scalarField(word("Qsum"), dict, p.size()))
      : std::move(scalarField(p.size(), 0.0))
    ),
    QsumNbr_
    (
        dict.found("QsumNbr")
      ? std::move(scalarField(word("QsumNbr"), dict, p.size()))
      : std::move(scalarField(p.size(), 0.0))
    ),
    qURF_(dict.lookupOrDefault<scalar>("Qurf", 0.3))
{
    if (dict.found("contactResistance") && dict.found("conductivityLayers"))
    {
        FatalErrorInFunction
            << "Either contact resistance or contact layers should be "
            << "specified in boundary "
            << p.name() << " of field " << iF.name() << endl;
    }

    if (dict.found("contactResistance"))
    {
        contactResistance_ = dict.lookup<scalar>("contactResistance");
        contactInsulance_ = contactResistance_*gSum(p.magSf());
    }
    else if (dict.found("conductivityLayers") || dict.found("thicknessLayers"))
    {
        conductivityLayers_ = dict.lookup<scalarList>("conductivityLayers");
        thicknessLayers_ = dict.lookup<scalarList>("thicknessLayers");

        if (thicknessLayers_.size() != conductivityLayers_.size())
        {
            FatalIOErrorInFunction(dict)
                << "conductivityLayers and thicknessLayers must be lists of "
                << "the same length" << nl << endl;
        }

        if (conductivityLayers_.size() > 0)
        {
            // Calculate effective contact insulance by harmonic averaging
            forAll(conductivityLayers_, iLayer)
            {
                contactInsulance_ +=
                    thicknessLayers_[iLayer] / conductivityLayers_[iLayer];
            }
        }
    }
    if (!restart())
    {
        // Prevent user initialising with a large initial heat flux to start
        scalarField::operator=(patchInternalField());
    }

    iF.mesh().schemes().setFluxRequired(iF.name());
}


regionCoupledFluxFvPatchScalarField::
regionCoupledFluxFvPatchScalarField
(
    const regionCoupledFluxFvPatchScalarField& ptf
)
:
    regionCoupledFvPatchField<scalar>(ptf),
    faceBCoeff_(ptf.faceBCoeff_),
    faceICoeff_(ptf.faceICoeff_),
    faceCorr_(ptf.faceCorr_),
    contactResistance_(ptf.contactResistance_),
    conductivityLayers_(ptf.conductivityLayers_),
    thicknessLayers_(ptf.thicknessLayers_),
    contactInsulance_(ptf.contactInsulance_),
    Qsum_(ptf.Qsum_),
    QsumNbr_(ptf.QsumNbr_),
    qURF_(ptf.qURF_)
{}


regionCoupledFluxFvPatchScalarField::
regionCoupledFluxFvPatchScalarField
(
    const regionCoupledFluxFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCoupledFvPatchField<scalar>(ptf, iF),
    faceBCoeff_(ptf.faceBCoeff_),
    faceICoeff_(ptf.faceICoeff_),
    faceCorr_(ptf.faceCorr_),
    contactResistance_(ptf.contactResistance_),
    conductivityLayers_(ptf.conductivityLayers_),
    thicknessLayers_(ptf.thicknessLayers_),
    contactInsulance_(ptf.contactInsulance_),
    Qsum_(ptf.Qsum_),
    QsumNbr_(ptf.QsumNbr_),
    qURF_(ptf.qURF_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void regionCoupledFluxFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    regionCoupledFvPatchField<scalar>::autoMap(m);
    m(faceBCoeff_, faceBCoeff_);
    m(faceICoeff_, faceICoeff_);
    m(faceCorr_, faceCorr_);
    m(Qsum_, Qsum_);
    m(QsumNbr_, QsumNbr_);
}


void regionCoupledFluxFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    regionCoupledFvPatchField<scalar>::rmap(ptf, addr);

    const regionCoupledFluxFvPatchScalarField& tiptf =
        refCast<const regionCoupledFluxFvPatchScalarField>(ptf);

    faceBCoeff_.rmap(tiptf.faceBCoeff_, addr);
    faceICoeff_.rmap(tiptf.faceICoeff_, addr);
    faceCorr_.rmap(tiptf.faceCorr_, addr);
    Qsum_.rmap(tiptf.Qsum_, addr);
    QsumNbr_.rmap(tiptf.QsumNbr_, addr);
}


void regionCoupledFluxFvPatchScalarField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    regionCoupledFvPatchField<scalar>::autoMapGIB(mapper);
    mapper.map(faceBCoeff_, scalar(0));
    mapper.map(faceICoeff_, scalar(0));
    mapper.map(faceCorr_, scalar(0));
    mapper.map(Qsum_, scalar(0));
    mapper.map(QsumNbr_, scalar(0));
}


tmp<scalarField> regionCoupledFluxFvPatchScalarField::snGrad() const
{
    if (!regionCoupled())
    {
        return regionCoupledFvPatchField<scalar>::snGrad();
    }

    return patch().deltaCoeffs() * (*this - patchInternalField());
}


tmp<scalarField>
regionCoupledFluxFvPatchScalarField::snGrad(const scalarField&) const
{
    return snGrad();
}


void regionCoupledFluxFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!regionCoupled())
    {
        regionCoupledFvPatchField<scalar>::evaluate();
        return;
    }

    updateCoeffs();
    forceAssign
    (
        faceICoeff_*patchInternalField()
      + faceBCoeff_*patchNeighbourField()
      + faceCorr_
    );
    fvPatchScalarField::evaluate();
}


void regionCoupledFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    fvPatchScalarField::updateCoeffs();
}


tmp<scalarField>
regionCoupledFluxFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    if (!regionCoupled())
    {
        return regionCoupledFvPatchField<scalar>::valueInternalCoeffs(w);
    }

    return tmp<scalarField>(new scalarField(patch().size(), 0));
}


tmp<scalarField>
regionCoupledFluxFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    if (!regionCoupled())
    {
        return regionCoupledFvPatchField<scalar>::valueBoundaryCoeffs(w);
    }

    // Note - for interfaces, this is a coefficient of the boundary value, not
    // the value itself
    return tmp<scalarField>(new scalarField(patch().size(), 1));
}


tmp<scalarField>
regionCoupledFluxFvPatchScalarField::gradientInternalCoeffs() const
{
    if (!regionCoupled())
    {
        return regionCoupledFvPatchField<scalar>::gradientInternalCoeffs();
    }

    return patch().deltaCoeffs()
      * (valueInternalCoeffs(patch().weights()) - 1.0);
}


tmp<scalarField>
regionCoupledFluxFvPatchScalarField::gradientBoundaryCoeffs() const
{
    if (!regionCoupled())
    {
        return regionCoupledFvPatchField<scalar>::gradientBoundaryCoeffs();
    }

    // Note - for interfaces, this is a coefficient of the boundary value,
    // not a source value
    return patch().deltaCoeffs() * valueBoundaryCoeffs(patch().weights());
}


void regionCoupledFluxFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    if (!regionCoupled())
    {
        return;
    }

    // Store matrix coefficients from here so they are available on both
    // sides when regionCoupledBoundaryCoeffs() is called
    const label patchi = this->patch().index();
    Cb_ = stabilise(fieldScale(matrix.boundaryCoeffs()[patchi]), VSMALL);
    Ci_ = matrix.internalCoeffs()[patchi];
    fluxCorr_ = scalarField(this->size(), Zero);
    if (matrix.faceFluxCorrectionPtr())
    {
        fluxCorr_ += matrix.faceFluxCorrectionPtr()->boundaryField()[patchi];
    }
}


void regionCoupledFluxFvPatchScalarField::regionCoupledBoundaryCoeffs
(
    const fvMatrix<scalar>& matrix,
    const scalarField& bouCoeffs,
    const scalarField& intCoeffs,
    scalarField& coupledBouCoeffs,
    scalarField& coupledIntCoeffs
)
{
    if (!regionCoupled())
    {
        return;
    }

    // Existing boundary/internal coeffs are the coefficients of the boundary
    // fluxes in the equation as a function of the (negative of the) face value
    // and internal cell value (bouCoeffs, intCoeffs). Multiply these by the
    // coefficients of the boundary values as a function of internal and
    // neighbour cell values (faceICoeff_, faceBCoeff_) to get coefficients of
    // boundary fluxes as a function of internal and neighbour cell values
    // (coupledBouCoeffs, coupledIntCoeffs).
    const regionCoupledFluxFvPatchScalarField& nbrPatchField =
        refCast<const regionCoupledFluxFvPatchScalarField>
        (
            neighbourFvPatchField()
        );
    const scalarField& magSf(patch().magSf());
    const scalarField& nbrMagSf(nbrPatchField.patch().magSf());

    scalarField CbNbrIns
    (
        scalar(1)
       /(
            scalar(1)
           /interpolateFromNeighbour
            (
                nbrPatchField.fluxScale(nbrPatchField.Cb_), fluxScale(Cb_)
            )
          + (contactInsulance_+nbrPatchField.contactInsulance_)/magSf
        )
    );

    scalarField weight(fluxScale(Cb_)/(fluxScale(Cb_) + CbNbrIns));

    faceICoeff_ = fieldScale(Ci_/Cb_*weight);
    coupledIntCoeffs = intCoeffs - bouCoeffs*faceICoeff_;

    faceBCoeff_ =
        fieldScale
        (
            interpolateFromNeighbour
            (
                nbrPatchField.Ci_/nbrPatchField.Cb_,
                Ci_/Cb_
            )*(1-weight)
        );
    coupledBouCoeffs = bouCoeffs*faceBCoeff_;

    // Update sources (on both sides), using Newton's method to account for
    // possible face value dependence

    label i = 0;
    scalar ppsi_err = GREAT, pNbrPsi_err = GREAT;
    scalar maxErr = 10*SMALL;
    label maxLoops = 100;

    scalarField ppsi_new(negFieldOffset(invFieldScale(*this)));
    scalarField pNbrPsi_new
    (
        interpolateFromNeighbour
        (
            nbrPatchField.negFieldOffset
            (
                nbrPatchField.invFieldScale(nbrPatchField)
            ),
            negFieldOffset(invFieldScale(*this))
        )
    );
    scalarField ppsi_old, pNbrPsi_old;

    scalarField Q, nbrQ;

    // We formulate everything such as to match flux with a fictional face
    // on the neighbour side which is coincident with the owner side.
    // Doing it exactly would require an NCC approach with sub-faces
    scalarField zeroField(size(), 0);
    scalarField nbrCi
    (
        interpolateFromNeighbour
        (
            nbrPatchField.fluxScale(nbrPatchField.Ci_)/nbrMagSf,
            // Default is zero to give zero flux
            zeroField
        )*magSf
    );
    scalarField nbrIntField
    (
        interpolateFromNeighbour
        (
            nbrPatchField.patchInternalField(),
            // Default doesn't matter since it will be multiplied by zero Ci
            zeroField
        )
    );
    scalarField nbrCb
    (
        interpolateFromNeighbour
        (
            nbrPatchField.fluxScale(nbrPatchField.Cb_)/nbrMagSf,
            // Default is zero to give zero flux
            zeroField
        )*magSf
    );
    scalarField nbrFluxCorr
    (
        interpolateFromNeighbour
        (
            nbrPatchField.fluxScale
            (
                nbrPatchField.fluxCorr_
            )/nbrMagSf,
            // Default is zero to give zero flux
            zeroField
        )*magSf
      - nbrCb
       *interpolateFromNeighbour
        (
            nbrPatchField.fieldOffset(scalarField(nbrPatchField.size(), 0)),
            zeroField
        )
    );

    do
    {
        ppsi_old = ppsi_new;
        pNbrPsi_old = pNbrPsi_new;
        scalarField bSourceDeriv, bnbrSourceDeriv;
        Q =
            fluxScale
            (
                boundarySources(fieldScale(fieldOffset(ppsi_old)), bSourceDeriv)
            );
        bSourceDeriv = fluxScale(fieldScale(bSourceDeriv));
        nbrQ =
            nbrPatchField.boundarySources
            (
                nbrPatchField.fieldScale
                (
                    nbrPatchField.fieldOffset
                    (
                        interpolateToNeighbour
                        (
                            pNbrPsi_old,
                            nbrPatchField.negFieldOffset
                            (
                                nbrPatchField.invFieldScale(nbrPatchField)
                            )
                        )
                    )
                ),
                bnbrSourceDeriv
            );
        nbrQ =
            interpolateFromNeighbour(nbrPatchField.fluxScale(nbrQ), zeroField);
        bnbrSourceDeriv =
            interpolateFromNeighbour
            (
                nbrPatchField.fluxScale
                (
                    nbrPatchField.fieldScale(bnbrSourceDeriv)
                ),
                zeroField
            );

        scalarField flux
        (
            fluxScale
            (
                Ci_*patchInternalField() - Cb_*ppsi_old + fluxCorr_
              - Cb_*fieldOffset(scalarField(size(), 0))
            )
        );
        scalarField nbrFlux
        (
            nbrCi*nbrIntField - nbrCb*pNbrPsi_old + nbrFluxCorr
        );

        // Outer flux balance equation
        scalarField f1(flux/magSf + Q + nbrFlux/magSf + nbrQ);
        scalarField df1dpsi(fluxScale(-Cb_/magSf)+bSourceDeriv);
        scalarField df1dPsiNbr(-nbrCb/magSf + bnbrSourceDeriv);

        // Layer correction equation
        scalarField f2
        (
            ppsi_old
          - contactInsulance_*(flux/magSf + Q)
          - pNbrPsi_old
          + nbrPatchField.contactInsulance_*(nbrFlux/magSf + nbrQ)
        );
        scalarField df2dpsi
        (
            scalarField(size(), scalar(1))
          - contactInsulance_*(fluxScale(-Cb_/magSf)+bSourceDeriv)
        );
        scalarField df2dpsiNbr
        (
           -scalarField(size(), scalar(1))
          + nbrPatchField.contactInsulance_*(-nbrCb/magSf + bnbrSourceDeriv)
        );

        // Solve (effectively block-Jacobi, since the interpolation to/from
        // the neighbouring side boundary source introduces some coupling
        // between neighbouring faces which we do not account for)

        scalarField D(df1dpsi*df2dpsiNbr - df1dPsiNbr*df2dpsi);
        ppsi_new =
            ppsi_old - (df2dpsiNbr*f1-df1dPsiNbr*f2)/stabilise(D, SMALL);
        pNbrPsi_new =
            pNbrPsi_old - (-df2dpsi*f1+df1dpsi*f2)/stabilise(D, SMALL);

        ppsi_err = gMax(mag(ppsi_new-ppsi_old)/stabilise(mag(ppsi_old), SMALL));
        pNbrPsi_err =
            gMax
            (
                mag(pNbrPsi_new-pNbrPsi_old)/stabilise(mag(pNbrPsi_old), SMALL)
            );

        if (debug)
        {
            Info<< "Patch " << this->patch().name() << ", Iter: " << i
                << ", Deltas: " << ppsi_err << ", " << pNbrPsi_err << endl;
        }

        i++;
    }
    while (i < maxLoops && (ppsi_err > maxErr || pNbrPsi_err > maxErr));

    if (i == maxLoops)
    {
        WarningInFunction
            << "Non-convergence in Newton's method for "
            << "temperature-dependent boundary source, patch: "
            << this->patch().name() << nl << endl;
    }

    // Apply under-relaxation to Qsum
    Qsum_ = qURF_*Q + (1-qURF_)*Qsum_;
    QsumNbr_ = qURF_*nbrQ + (1-qURF_)*QsumNbr_;

    // Compute final boundary face value correction
    // Keeping copy so that correct value can be reproduced in evaluate()
    faceCorr_ =
        fieldScale
        (
            (fluxCorr_ + Qsum_*magSf)/Cb_*weight
          + fieldOffset((nbrFluxCorr + QsumNbr_*magSf)/stabilise(nbrCb, VSMALL))*(1-weight)
        );
}


void regionCoupledFluxFvPatchScalarField::write(Ostream& os) const
{
    regionCoupledFvPatchField<scalar>::write(os);
    if (conductivityLayers_.size())
    {
        conductivityLayers_.writeEntry("conductivityLayers", os);
        thicknessLayers_.writeEntry("thicknessLayers", os);
    }
    else
    {
        os.writeEntry("contactResistance", contactResistance_);
    }

    faceBCoeff_.writeEntry("faceBCoeff", os);
    faceICoeff_.writeEntry("faceICoeff", os);
    faceCorr_.writeEntry("faceCorr", os);
    Qsum_.writeEntry("Qsum", os);
    QsumNbr_.writeEntry("QsumNbr", os);

    writeEntryIfDifferent<scalar>
    (
        os,
        "Qurf",
        0.3,
        qURF_
    );

    this->writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


makePatchTypeField
(
    fvPatchScalarField,
    regionCoupledFluxFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
