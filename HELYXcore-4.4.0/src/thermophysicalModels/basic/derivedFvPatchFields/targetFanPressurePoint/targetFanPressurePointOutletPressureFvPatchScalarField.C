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
    (c) 2020-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "derivedFvPatchFields/targetFanPressurePoint/targetFanPressurePointOutletPressureFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "cfdTools/general/solutionControl/helyxCoupledControl/helyxCoupledControl.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::targetFanPressurePointOutletPressureFvPatchScalarField::checkPhiDimensions
(
    const fvsPatchField<scalar>& phiP
) const
{
    const dimensionSet& phiDims = phiP.internalField().dimensions();
    const bool isPhiMassFlow = (phiDims == dimDensity*dimVelocity*dimArea);
    if (!isPhiMassFlow && (phiDims != dimVelocity*dimArea))
    {
        FatalErrorInFunction
            << phiName_ << " field dimensions are not consistent"
            << abort(FatalError);
    }
    return isPhiMassFlow;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::targetFanPressurePointOutletPressureFvPatchScalarField::
targetFanPressurePointOutletPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    refPoint_(vector::zero),
    initialP_(0),
    isTimeCurve_(false),
    fanCurve_(new Function1Types::Constant<scalar>("fanCurve", 0.0)),
    timeFanCurve_(nullptr),
    isInitial_(true)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::targetFanPressurePointOutletPressureFvPatchScalarField::
targetFanPressurePointOutletPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    refPoint_(dict.lookup<point>("referencePoint")),
    initialP_(dict.lookupOrDefault<scalar>("initialPressure", 0)),
    isInitial_(true)
{
    frameName_ = dict.lookupOrDefault<word>("referenceFrame", "");
    inletFlux_ = true;

    isTimeCurve_ = dict.lookupOrDefault("timeFanCurve", false);
    if (isTimeCurve_)
    {
        timeFanCurve_ = Function2<scalar>::New("fanCurve", dict);
    }
    else
    {
        if (dict.found("fanCurve"))
        {
            fanCurve_ = Function1<scalar>::New("fanCurve", dict);
        }
        else
        {
            FatalIOErrorIn
            (
                "fanPressureFvPatchScalarField::"
                "fanPressureFvPatchScalarField"
                "(const fvPatch&, const DimensionedField<vector, volMesh>&,"
                " const dictionary&)",
                dict
            )   << "Please supply 'fanCurve'"
                << exit(FatalIOError);
        }
    }

    if (dict.found("value"))
    {
        this->refValue() = scalarField("value", dict, p.size());
        fvPatchField<scalar>::operator= (this->refValue());
    }
    else
    {
        this->refValue() = 0.0;
    }
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::targetFanPressurePointOutletPressureFvPatchScalarField::
targetFanPressurePointOutletPressureFvPatchScalarField
(
    const targetFanPressurePointOutletPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    refPoint_(ptf.refPoint_),
    initialP_(ptf.initialP_),
    isTimeCurve_(ptf.isTimeCurve_),
    fanCurve_(!isTimeCurve_ ? ptf.fanCurve_().clone().ptr() : nullptr),
    timeFanCurve_(isTimeCurve_ ? ptf.timeFanCurve_().clone().ptr() : nullptr),
    isInitial_(ptf.isInitial_)
{}


Foam::targetFanPressurePointOutletPressureFvPatchScalarField::
targetFanPressurePointOutletPressureFvPatchScalarField
(
    const targetFanPressurePointOutletPressureFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    refPoint_(ptf.refPoint_),
    initialP_(ptf.initialP_),
    isTimeCurve_(ptf.isTimeCurve_),
    fanCurve_(!isTimeCurve_ ? ptf.fanCurve_().clone().ptr() : nullptr),
    timeFanCurve_(isTimeCurve_ ? ptf.timeFanCurve_().clone().ptr() : nullptr),
    isInitial_(ptf.isInitial_)
{}


Foam::targetFanPressurePointOutletPressureFvPatchScalarField::
targetFanPressurePointOutletPressureFvPatchScalarField
(
    const targetFanPressurePointOutletPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    refPoint_(ptf.refPoint_),
    initialP_(ptf.initialP_),
    isTimeCurve_(ptf.isTimeCurve_),
    fanCurve_(!isTimeCurve_ ? ptf.fanCurve_().clone().ptr() : nullptr),
    timeFanCurve_(isTimeCurve_ ? ptf.timeFanCurve_().clone().ptr() : nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::targetFanPressurePointOutletPressureFvPatchScalarField::referenceValue
(
    const word& fieldName
)
{
    const volScalarField& ptf = db().lookupObject<volScalarField>(fieldName);
    label cellI =  ptf.mesh().findNearestCell(refPoint_);

    List<scalar> cellDist(Pstream::nProcs());
    cellDist[Pstream::myProcNo()] = GREAT;
    List<scalar> refValue(Pstream::nProcs());
    refValue[Pstream::myProcNo()] = -1;
    scalar& myCellDist = cellDist[Pstream::myProcNo()];

    if (cellI>-1)
    {
        myCellDist = mag(ptf.mesh().cellCentres()[cellI] -refPoint_);
        refValue[Pstream::myProcNo()] = ptf.internalField()[cellI];
    }

    Pstream::allGatherList(cellDist);
    Pstream::allGatherList(refValue);
    scalar closestCellDist = GREAT;
    scalar closestRefValue(0);
    forAll(cellDist, proci)
    {
        if (cellDist[proci] < closestCellDist)
        {
            closestCellDist = cellDist[proci];
            closestRefValue = refValue[proci];
        }
    }
    if (debug)
    {
        Pout<< "Found closest cell at distance " << closestCellDist
             << ", ref. value " << closestRefValue
             << ", for field " << fieldName << endl;
    }

    return closestRefValue;
}


void Foam::targetFanPressurePointOutletPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (this->internalField().name() == "pPot")
    {
        this->valueFraction() = 1.0;
        this->refGrad() = 0.0;
    }
    else
    {
       label numIters = 20;
       for (label iter=0; iter<numIters; iter++)
       {

            fvsPatchField <scalar> phiP =
                patch().lookupPatchFieldInDb<surfaceScalarField, scalar>
                (
                    this->db(),
                    phiName_
                );
            const bool isPhiMassFlow = checkPhiDimensions(phiP);

            // Load rho for compressible
            const fvPatchField <scalar> *rhoP =
                isPhiMassFlow
                ? &patch().lookupPatchFieldInDb<volScalarField, scalar>
                (
                    this->db(),
                    this->rhoName_
                )
                : nullptr;

            if (!this->db().foundObject<helyxCoupledControl>(solutionControl::typeName)) {
                if (isPhiMassFlow)
                {
                    phiP = phiP/(*rhoP);
                    makeRelative(phiP);
                    // scalarField relPhiP(phiP.size(), 0.0);
                    // makeRelative(relPhiP);
                    // phiP += (*rhoP) * relPhiP;
                }
                else
                {
                    makeRelative(phiP);
                }
            }

            const scalar currentFlux = gSum(phiP);

            scalar pdFan;
            if (timeFanCurve_)
            {
                pdFan =
                    timeFanCurve_->value
                    (
                        this->db().time().timeOutputValue(),
                        max(currentFlux, 0.0)
                    );
            }
            else
            {
                pdFan = fanCurve_->value(max(currentFlux, 0.0));
            }

            scalarField &p = *this;

            scalar pAvg
            (
                gSum((*this)*(*this).patch().magSf())/gSum(patch().magSf())
            );

            scalar referencePressure = referenceValue("p");

            scalar dp = -(pdFan-mag(referencePressure));

            if (isInitial_)
            {
                forAll(p, fI)
                {
                    p[fI] = initialP_;
                }
            }
            isInitial_ = false;
            p += dp;

            if (iter == numIters-1)
            {
                Info<< "Current flux: " << currentFlux
                     << ", p from table: " << pdFan
                     << ", p target: " << referencePressure
                     << ", dp: " << dp
                     << ", initial p (outlet): " << pAvg
                     << ", new p (outlet): " << gAverage(p) << endl;
            }

            this->refValue() = p;
            this->valueFraction() = 1.0;

       } //iter

    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::targetFanPressurePointOutletPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);

    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    os.writeEntry("referencePoint", refPoint_);
    os.writeEntry("initialPressure", initialP_);
    os.writeEntry("timeFanCurve", isTimeCurve_);
    if (isTimeCurve_)
    {
        timeFanCurve_->writeData(os);
    }
    else
    {
        fanCurve_->writeData(os);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        targetFanPressurePointOutletPressureFvPatchScalarField
    );
}

// ************************************************************************* //
