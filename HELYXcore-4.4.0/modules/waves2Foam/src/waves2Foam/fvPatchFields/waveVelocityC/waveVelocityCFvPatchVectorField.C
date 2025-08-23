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
    (c) held by original author

\*---------------------------------------------------------------------------*/

#include "waveVelocityCFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#if EXTBRANCH==1
    #if 310<OFVERSION
        #include "foamTime.H"
    #else
        #include "db/Time/Time.H"
    #endif
#elif OFPLUSBRANCH==1
    #include "db/Time/Time.H"
#else
    #include "db/Time/Time.H"
#endif

#include "dimensionedTypes/dimensionedVector/dimensionedVector.H"
#include "volMesh/volMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

waveVelocityCFvPatchVectorField::waveVelocityCFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchField<vector>(p, iF),
    #if EXTBRANCH==1
        convexPolyhedral(this->dimensionedInternalField().mesh(), true),
    #elif OFPLUSBRANCH==1
        #if OFVERSION<1706
            convexPolyhedral(this->dimensionedInternalField().mesh(), true),
        #else
            convexPolyhedral(this->internalField().mesh(), true),
        #endif
    #else
        #if OFVERSION<400
            convexPolyhedral(this->dimensionedInternalField().mesh(), true),
        #else
            convexPolyhedral(this->internalField().mesh(), true),
        #endif
    #endif
    waveProps_(waveTheories::waveTheory::New(this->patch().name(),
    #if EXTBRANCH==1
                this->dimensionedInternalField().mesh()
    #elif OFPLUSBRANCH==1
        #if OFVERSION<1706
                this->dimensionedInternalField().mesh()
        #else
                this->internalField().mesh()
        #endif
    #else
        #if OFVERSION<400
                this->dimensionedInternalField().mesh()
        #else
                this->internalField().mesh()
        #endif
    #endif
    )),
    phiName_("phi"),
    thresholdValue(),
    thresholdAxis()
{
    this->refValue() = pTraits<vector>::zero;
    this->refGrad() = pTraits<vector>::zero;
    this->valueFraction() = 0.0;

}


waveVelocityCFvPatchVectorField::waveVelocityCFvPatchVectorField
(
    const waveVelocityCFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<vector>(ptf, p, iF, mapper),
    #if EXTBRANCH==1
        convexPolyhedral(this->dimensionedInternalField().mesh(), true),
    #elif OFPLUSBRANCH==1
        #if OFVERSION<1706
            convexPolyhedral(this->dimensionedInternalField().mesh(), true),
        #else
            convexPolyhedral(this->internalField().mesh(), true),
        #endif
    #else
        #if OFVERSION<400
            convexPolyhedral(this->dimensionedInternalField().mesh(), true),
        #else
            convexPolyhedral(this->internalField().mesh(), true),
        #endif
    #endif
    waveProps_(ptf.waveProps_),
    phiName_(ptf.phiName_),
    thresholdValue(ptf.thresholdValue),
    thresholdAxis(ptf.thresholdAxis)

{}


waveVelocityCFvPatchVectorField::waveVelocityCFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<vector>(p, iF, dict, false),
    #if EXTBRANCH==1
        convexPolyhedral(this->dimensionedInternalField().mesh(), true),
    #elif OFPLUSBRANCH==1
        #if OFVERSION<1706
            convexPolyhedral(this->dimensionedInternalField().mesh(), true),
        #else
            convexPolyhedral(this->internalField().mesh(), true),
        #endif
    #else
        #if OFVERSION<400
            convexPolyhedral(this->dimensionedInternalField().mesh(), true),
        #else
            convexPolyhedral(this->internalField().mesh(), true),
        #endif
    #endif
    waveProps_(waveTheories::waveTheory::New(this->patch().name(),
    #if EXTBRANCH==1
                this->dimensionedInternalField().mesh()
    #elif OFPLUSBRANCH==1
        #if OFVERSION<1706
                this->dimensionedInternalField().mesh()
        #else
                this->internalField().mesh()
        #endif
    #else
        #if OFVERSION<400
                this->dimensionedInternalField().mesh()
        #else
                this->internalField().mesh()
        #endif
    #endif
    )),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    thresholdValue(dict.lookup<scalar>("thresholdValue")),
    thresholdAxis(dict.lookup("thresholdAxis"))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = Zero;
}


waveVelocityCFvPatchVectorField::waveVelocityCFvPatchVectorField
(
    const waveVelocityCFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchField<vector>(ptf, iF),
    #if EXTBRANCH==1
        convexPolyhedral(this->dimensionedInternalField().mesh(), true),
    #elif OFPLUSBRANCH==1
        #if OFVERSION<1706
            convexPolyhedral(this->dimensionedInternalField().mesh(), true),
        #else
            convexPolyhedral(this->internalField().mesh(), true),
        #endif
    #else
        #if OFVERSION<400
            convexPolyhedral(this->dimensionedInternalField().mesh(), true),
        #else
            convexPolyhedral(this->internalField().mesh(), true),
        #endif
    #endif
    waveProps_(waveTheories::waveTheory::New(this->patch().name(),
    #if EXTBRANCH==1
                this->dimensionedInternalField().mesh()
    #elif OFPLUSBRANCH==1
        #if OFVERSION<1706
                this->dimensionedInternalField().mesh()
        #else
                this->internalField().mesh()
        #endif
    #else
        #if OFVERSION<400
                this->dimensionedInternalField().mesh()
        #else
                this->internalField().mesh()
        #endif
    #endif
    )),
    phiName_(ptf.phiName_),
    thresholdValue(ptf.thresholdValue),
    thresholdAxis(ptf.thresholdAxis)

{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void waveVelocityCFvPatchVectorField::signedPointToSurfaceDistance
(
    const pointField& pp,
    scalarField& sd
)
{
    forAll(pp, pointi)
    {
        sd[pointi] = signedPointToSurfaceDistance(pp[pointi]);
    }
}


scalar waveVelocityCFvPatchVectorField::signedPointToSurfaceDistance
(
    const point& pp
) const
{
    scalar temp = waveProps_->eta(pp, db().time().value() );
    temp += ( waveProps_->returnDir() & pp );
    temp *= -1.0;

    return temp;
}


// Update the coefficients associated with the patch field
void waveVelocityCFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    #if EXTBRANCH==1
        const fvMesh& mesh = this->dimensionedInternalField().mesh();
    #elif OFPLUSBRANCH==1
        #if OFVERSION<1706
            const fvMesh& mesh = this->dimensionedInternalField().mesh();
        #else
            const fvMesh& mesh = this->internalField().mesh();
        #endif
    #else
        #if OFVERSION<400
            const fvMesh& mesh = this->dimensionedInternalField().mesh();
        #else
            const fvMesh& mesh = this->internalField().mesh();
        #endif
    #endif
    const word patchName = this->patch().name();
    const label patchID = mesh.boundaryMesh().findPatchID(patchName);
    const scalarField& magSf( mesh.magSf().boundaryField()[patchID] );

    const label start = patch().patch().start();

    // Loop through all boundary faces on the current patch and overwrite the boundary conditions
    forAll(magSf, facei)
    {

        localFace lf = this->divideFace(facei + start);

        // If the face is below the instantaneous wave elevation,
        // change face to fixedValue based on wave velocity
        if (lf.isNegFace())
        {
            this->refValue()[facei]
                = waveProps_->U( lf.negCentre(), db().time().value() );
            this->valueFraction()[facei] = 1.0;
        }
        else
        {
            // Save the axis as an integer index
            int ax(2);
            if (thresholdAxis == "x")
            {
                ax = 0;
            }
            else if (thresholdAxis == "y")
            {
                ax = 1;
            }
            else if (thresholdAxis == "z")
            {
                ax = 2;
            }
            else
            {
                FatalErrorInFunction
                    << "Incorrect value entered for thresholdAxis" << thresholdAxis
                    << "Correct values are x, y, and z"
                    << exit(FatalError);
            }

            // If the face is above the instantaneous wave elevation but below the specified threshold,
            // change face to fixedValue based on wind velocity
            if (lf.posCentre()[ax] < thresholdValue)
            {
                this->refValue()[facei]
                    = waveProps_->windVelocity( db().time().value() );
                this->valueFraction()[facei] = 1.0;
            }

            // If the face is above the instantaneous wave elevation and above the specified threshold,
            // change face to zeroGradient
            else
            {
                this->refGrad()[facei]       = vector::zero;
                this->refValue()[facei] = waveProps_->windVelocity( db().time().value() );
                this->valueFraction()[facei] = 0.0;
            }
        }
    }

    mixedFvPatchField<vector>::updateCoeffs();
    evaluate();
}


// Evaluate the field on the patch
void waveVelocityCFvPatchVectorField::evaluate()
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<vector>::operator=
    (
        this->valueFraction()*this->refValue()
      +
        (1.0 - this->valueFraction())*
        (
            this->patchInternalField()
          + this->refGrad()/this->patch().deltaCoeffs()
        )
    );

    fvPatchField<vector>::evaluate();
}


// Write
void waveVelocityCFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    this->refValue().writeEntry("refValue", os);
    this->refGrad().writeEntry("refGradient", os);
    this->valueFraction().writeEntry("valueFraction", os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    os.writeEntry<scalar>("thresholdValue", thresholdValue);
    os.writeEntry<word>("thresholdAxis", thresholdAxis);

    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, waveVelocityCFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
