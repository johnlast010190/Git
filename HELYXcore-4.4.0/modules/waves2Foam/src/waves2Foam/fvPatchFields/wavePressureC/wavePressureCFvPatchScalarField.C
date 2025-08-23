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
    (c) 2011-2018 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "wavePressureCFvPatchScalarField.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/surfaceFields/surfaceFields.H"

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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wavePressureCFvPatchScalarField::
wavePressureCFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
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
    p0_(p.size(), 0.0),
    thresholdValue(),
    thresholdAxis()
{}


Foam::wavePressureCFvPatchScalarField::
wavePressureCFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
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
    p0_("p0", dict, p.size()),
    thresholdValue(dict.lookup<scalar>("thresholdValue")),
    thresholdAxis(dict.lookup("thresholdAxis"))
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }

    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 0.0;
}


Foam::wavePressureCFvPatchScalarField::
wavePressureCFvPatchScalarField
(
    const wavePressureCFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
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
    p0_(mapper(ptf.p0_)),
    thresholdValue(ptf.thresholdValue),
    thresholdAxis(ptf.thresholdAxis)
{}


Foam::wavePressureCFvPatchScalarField::
wavePressureCFvPatchScalarField
(
    const wavePressureCFvPatchScalarField& wbppsf
)
:
    mixedFvPatchScalarField(wbppsf),
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
    p0_(wbppsf.p0_),
    thresholdValue(wbppsf.thresholdValue),
    thresholdAxis(wbppsf.thresholdAxis)
{}


Foam::wavePressureCFvPatchScalarField::
wavePressureCFvPatchScalarField
(
    const wavePressureCFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wbppsf, iF),
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
    p0_(wbppsf.p0_),
    thresholdValue(wbppsf.thresholdValue),
    thresholdAxis(wbppsf.thresholdAxis)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void wavePressureCFvPatchScalarField::signedPointToSurfaceDistance
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

scalar wavePressureCFvPatchScalarField::signedPointToSurfaceDistance
(
    const point& pp
) const
{
    scalar temp = waveProps_->eta(pp, db().time().value() );
    temp += ( waveProps_->returnDir() & pp );
    temp *= -1.0;

    return temp;
}

void Foam::wavePressureCFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Read in the mesh
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

    // Initially set everything to fixedValue boundary condition
    valueFraction() = 1;
    refValue() = p0_;
    refGrad() = Zero;

    // Loop through all boundary faces on the current patch and overwrite the boundary conditions
    forAll(magSf, facei)
    {

        localFace lf = this->divideFace(facei + start);

        // If the face is below the instantaneous wave elevation,
        // change face to zeroGradient
        if (lf.isNegFace())
        {
            this->valueFraction()[facei] = 0.0;
            this->refValue()[facei] = 0.0;
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
            // change face to zeroGradient
            if (lf.posCentre()[ax] < thresholdValue)
            {
                this->valueFraction()[facei] = 0.0;
                this->refValue()[facei] = 0.0;
            }
        }
    }

    mixedFvPatchField<scalar>::updateCoeffs();
    mixedFvPatchField<scalar>::evaluate();
}

void Foam::wavePressureCFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry<scalar>("thresholdValue", thresholdValue);
    os.writeEntry<word>("thresholdAxis", thresholdAxis);
    p0_.writeEntry("p0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        wavePressureCFvPatchScalarField
    );
}

// ************************************************************************* //
