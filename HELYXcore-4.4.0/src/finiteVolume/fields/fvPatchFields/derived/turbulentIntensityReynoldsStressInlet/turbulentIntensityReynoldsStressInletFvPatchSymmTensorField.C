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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2010-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/turbulentIntensityReynoldsStressInlet/turbulentIntensityReynoldsStressInletFvPatchSymmTensorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentIntensityReynoldsStressInletFvPatchSymmTensorField::
turbulentIntensityReynoldsStressInletFvPatchSymmTensorField
(
    const fvPatch& p,
    const DimensionedField<symmTensor, volMesh>& iF
)
:
    inletOutletFvPatchSymmTensorField(p, iF),
    intensity_(0.0),
    UName_(IOobject::groupName("U", internalField().group())),
    rhoName_("rho")
{
    this->refValue() = pTraits<symmTensor>::zero;
    this->refGrad()  = pTraits<symmTensor>::zero;
    this->valueFraction() = 0.0;
}

Foam::turbulentIntensityReynoldsStressInletFvPatchSymmTensorField::
turbulentIntensityReynoldsStressInletFvPatchSymmTensorField
(
    const turbulentIntensityReynoldsStressInletFvPatchSymmTensorField& ptf,
    const fvPatch& p,
    const DimensionedField<symmTensor, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchSymmTensorField(ptf, p, iF, mapper),
    intensity_(ptf.intensity_),
    UName_(ptf.UName_),
    rhoName_(ptf.rhoName_)
{
}

Foam::turbulentIntensityReynoldsStressInletFvPatchSymmTensorField::
turbulentIntensityReynoldsStressInletFvPatchSymmTensorField
(
    const fvPatch& p,
    const DimensionedField<symmTensor, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchSymmTensorField(p, iF),
    intensity_(dict.lookup<scalar>("intensity")),
    UName_
    (
        dict.lookupOrDefault<word>
        (
            "U", IOobject::groupName("U", internalField().group())
        )
    ),
    rhoName_
    (
        dict.lookupOrDefault<word>
        (
            "rho",
            internalField().group() == word::null
          ? "rho"
          : IOobject::groupName("thermo:rho", internalField().group())
        )
    )
{

    this->phiName_ = dict.lookupOrDefault<word>("phi", word::null);

    if (intensity_ < 0 || intensity_ > 1)
    {
    FatalErrorInFunction
        << "Turbulence intensity should be specified as a fraction 0-1 "
        "of the mean velocity\n"
        "    value given is " << intensity_ << nl
        << "    on patch " << this->patch().name()
        << " of field " << this->internalField().name()
        << " in file " << this->internalField().objectPath()
        << exit(FatalError);
    }

    fvPatchSymmTensorField::operator=(symmTensorField("value", dict, p.size()));

    this->refValue() = *this;
    this->refGrad() = pTraits<symmTensor>::zero;
    this->valueFraction() = SMALL;
}

Foam::turbulentIntensityReynoldsStressInletFvPatchSymmTensorField::
turbulentIntensityReynoldsStressInletFvPatchSymmTensorField
(
    const turbulentIntensityReynoldsStressInletFvPatchSymmTensorField& ptf
)
:
    inletOutletFvPatchSymmTensorField(ptf),
    intensity_(ptf.intensity_),
    UName_(ptf.UName_),
    rhoName_(ptf.rhoName_)
{
}


Foam::turbulentIntensityReynoldsStressInletFvPatchSymmTensorField::
turbulentIntensityReynoldsStressInletFvPatchSymmTensorField
(
    const turbulentIntensityReynoldsStressInletFvPatchSymmTensorField& ptf,
    const DimensionedField<symmTensor, volMesh>& iF
)
:
    inletOutletFvPatchSymmTensorField(ptf, iF),
    intensity_(ptf.intensity_),
    UName_(ptf.UName_),
    rhoName_(ptf.rhoName_)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::turbulentIntensityReynoldsStressInletFvPatchSymmTensorField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalarField UpSqr
    (
        magSqr
        (
            patch().lookupPatchFieldInDb<volVectorField, vector>(db(), UName_)
        )
    );

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(this->detectPhiName(phiName_));

    const fvsPatchScalarField& phip =
        patch().lookupPatchFieldInDb<surfaceScalarField, scalar>
        (
            db(), this->detectPhiName(this->phiName_)
        );

    if (phi.dimensions() == dimArea*dimVelocity)
    {
        UpSqr = max(UpSqr, magSqr(phip/patch().magSf()));
    }
    else if (phi.dimensions() == dimDensity*dimArea*dimVelocity)
    {
        if (db().foundObject<volScalarField>(rhoName_))
        {
            const fvPatchScalarField& rhop =
                patch().lookupPatchFieldInDb<volScalarField, scalar>
                (
                    db(), rhoName_
                );

            UpSqr = max(UpSqr, magSqr(phip/rhop/patch().magSf()));
        }
        else
        {
            WarningInFunction
                << "Could not find density field " << rhoName_
                << " needed to calculate face flux velocity."
                << " Reverting to default behaviour." << endl;
        }
    }
    else
    {
        //default behaviour
    }

    //k=3/2*I^2*U^2
    //Rij=2/3*k*dij
    scalarField Rii( max(SMALL, (sqr(intensity_)*UpSqr)) );
    this->refValue() = Rii*symmTensor(1,0,0,1,0,1);
    this->valueFraction() = 1.0 - pos0(phip);

    inletOutletFvPatchSymmTensorField::updateCoeffs();
}


void Foam::turbulentIntensityReynoldsStressInletFvPatchSymmTensorField::write
(
    Ostream& os
) const
{
    fvPatchSymmTensorField::write(os);
    os.writeEntry("intensity", intensity_);
    writeEntryIfDifferent<word>
    (
        os, "U", IOobject::groupName("U", internalField().group()), UName_
    );
    writeEntryIfDifferent<word>(os, "phi", word::null, this->phiName_);
    writeEntryIfDifferent<word>
    (
        os,
        "rho",
        internalField().group() == word::null
      ? "rho"
      : IOobject::groupName("thermo:rho", internalField().group()),
        rhoName_
    );
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::turbulentIntensityReynoldsStressInletFvPatchSymmTensorField::
operator=
(
    const fvPatchField<symmTensor>& ptf
)
{
    fvPatchField<symmTensor>::operator=
    (
        this->valueFraction()*this->refValue()
        + (1 - this->valueFraction())*ptf
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchSymmTensorField,
        turbulentIntensityReynoldsStressInletFvPatchSymmTensorField
    );
}

// ************************************************************************* //
