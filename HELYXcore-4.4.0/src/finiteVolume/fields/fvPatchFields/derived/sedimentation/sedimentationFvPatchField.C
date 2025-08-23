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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "sedimentationFvPatchField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "global/unitConversion/unitConversion.H"
#include "cfdTools/general/include/fvCFD.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makeRemovablePatchTypeField
    (
        fvPatchScalarField,
        sedimentationFvPatchField
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sedimentationFvPatchField::sedimentationFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    constrainDisplacement_(false),
    sedimentationSurfaceConstrainName_(word::null),
    sedimentationSurfaceConstrain_(nullptr),
    normal_(0, 0, 1),
    timeIndex_(-1),
    minValue_(0.0),
    maxValue_(GREAT)
{}


Foam::sedimentationFvPatchField::sedimentationFvPatchField
(
    const sedimentationFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    constrainDisplacement_(ptf.constrainDisplacement_),
    sedimentationSurfaceConstrainName_(ptf.sedimentationSurfaceConstrainName_),
    sedimentationSurfaceConstrain_(ptf.sedimentationSurfaceConstrain_),
    normal_(ptf.normal_),
    timeIndex_(ptf.timeIndex_),
    minValue_(ptf.minValue_),
    maxValue_(ptf.maxValue_)
{}


Foam::sedimentationFvPatchField::sedimentationFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    constrainDisplacement_
    (
        dict.lookupOrDefault<Switch>("constrainDisplacement", false)
    ),
    sedimentationSurfaceConstrainName_
    (
        constrainDisplacement_
      ? dict.lookup<word>("surfaceConstrainName")
      : word::null
    ),
    sedimentationSurfaceConstrain_(nullptr),
    normal_(0, 0, 1),
    timeIndex_(-1),
    minValue_(dict.lookupOrDefault<scalar>("minValue", 0.0)),
    maxValue_(dict.lookupOrDefault<scalar>("maxValue", GREAT))
{
    const Time& runTime = this->patch().boundaryMesh().mesh().time();
    const objectRegistry& obr = this->db();

    // Load normal direction (gravity vector is prefered)
    if (dict.found("upDirection"))
    {
        normal_ = -dict.lookup<vector>("upDirection");
    }
    else if (obr.foundObject<uniformDimensionedVectorField>("g"))
    {
        normal_ =
            -obr.lookupObject<uniformDimensionedVectorField>("g").value();
    }
    else if
    (
        IOobject
        (
            "g",
            runTime.constant(),
            obr,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ).headerOk()
    )
    {
        normal_ =
            -uniformDimensionedVectorField
            (
                IOobject
                (
                    "g",
                    runTime.constant(),
                    obr,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            ).value();
    }
    else
    {
        WarningInFunction
            << "Boundary "
            << sedimentationFvPatchField::typeName
            << " didn't find gravity vector g or upDirection"
            << " in the dictionary. Using default value (0, 0, 1)."
            << nl;

    }
    normal_ = normalised(normal_);
}


Foam::sedimentationFvPatchField::sedimentationFvPatchField
(
    const sedimentationFvPatchField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf),
    constrainDisplacement_(ptf.constrainDisplacement_),
    sedimentationSurfaceConstrainName_(ptf.sedimentationSurfaceConstrainName_),
    sedimentationSurfaceConstrain_(ptf.sedimentationSurfaceConstrain_),
    normal_(ptf.normal_),
    timeIndex_(ptf.timeIndex_),
    minValue_(ptf.minValue_),
    maxValue_(ptf.maxValue_)
{}


Foam::sedimentationFvPatchField::sedimentationFvPatchField
(
    const sedimentationFvPatchField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
    constrainDisplacement_(ptf.constrainDisplacement_),
    sedimentationSurfaceConstrainName_(ptf.sedimentationSurfaceConstrainName_),
    sedimentationSurfaceConstrain_(ptf.sedimentationSurfaceConstrain_),
    normal_(ptf.normal_),
    timeIndex_(ptf.timeIndex_),
    minValue_(ptf.minValue_),
    maxValue_(ptf.maxValue_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sedimentationFvPatchField::~sedimentationFvPatchField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sedimentationFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const objectRegistry& obr = this->db();
    const Time& runTime = obr.time();

    scalarField values(*this);

    // Turbulent viscosity
    scalarField nut(patch().lookupPatchField<volScalarField, scalar>("nut"));

    // Set min nut to avoid division by zero
    forAll(nut, facei)
    {
        nut[facei] = max(nut[facei], 1e-8);
    }

    // Mass exchange field
    const fvPatchField<scalar>& Mp =
        patch().lookupPatchField<volScalarField, scalar>("M");

    // Limiting sedimentation by the surface beyond which the mesh can't move
    // anymore
    boolList belowTrashold(patch().size(), false);

    if (constrainDisplacement_)
    {
        const pointField& Cfs = patch().Cf();
        // Load the surface first time it is needed
        if (!sedimentationSurfaceConstrain_.valid())
        {
            sedimentationSurfaceConstrain_.reset
            (
                new triSurfaceMesh
                (
                    IOobject
                    (
                        sedimentationSurfaceConstrainName_,
                        runTime.constant(),
                        "triSurface",
                        obr,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    )
                )
            );
        }

        // Find all intersections with the surface
        vectorField nearestPoints
        (
            sedimentationSurfaceConstrain_().findRayNearest(Cfs, -normal_)
        );

        forAll(*this, facei)
        {
            belowTrashold[facei] =
                nearestPoints[facei] == vector(GREAT, GREAT, GREAT)
             || (Cfs[facei] & normal_) < (nearestPoints[facei] & normal_);
        }
    }

    const label curTimeIndex =
        runTime.subCycling()
      ? runTime.prevTimeState().timeIndex()
      : runTime.timeIndex();

    scalarField totalSedimentDepth
    (
        patch().lookupPatchField<volScalarField, scalar>("totalSedimentDepth")
    );

    if (curTimeIndex != timeIndex_)
    {
        // Update time index
        timeIndex_ = curTimeIndex;

        // Get the sedimentation drift at the cell center
        // adjacent to the boundary surface.
        const tmp<scalarField> internalValues((*this).patchInternalField());
        forAll(*this, facei)
        {
            // Concentration gradient in the z-direction
            const scalar divisor = (-patch().Sf()[facei] & normal_)*nut[facei];

            scalar gradient =
                mag(divisor) <= SMALL ? 0.0 : Mp[facei]/divisor;
            if
            (
                gradient > 0
             || belowTrashold[facei]
             || totalSedimentDepth[facei] <= VSMALL
            )
            {
                // Gradient is 0 when deposition occurs.
                gradient = 0.0;
            }

            values[facei] =
                internalValues()[facei]
              - gradient/patch().deltaCoeffs()[facei];
            values[facei] = max(min(values[facei], maxValue_), minValue_);
        }
    }

    this->forceAssign(values);

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void Foam::sedimentationFvPatchField::write(Ostream& os) const
{
    fixedValueFvPatchField<scalar>::write(os);

    os.writeEntry("constrainDisplacement", constrainDisplacement_);
    os.writeEntryIfDifferent<word>
    (
        "surfaceConstrainName",
        word::null,
        sedimentationSurfaceConstrainName_
    );
    os.writeEntryIfDifferent<vector>("upDirection", vector(0, 0, 1), normal_);
    os.writeEntryIfDifferent<scalar>("minValue", 0.0, minValue_);
    os.writeEntryIfDifferent<scalar>("maxValue", GREAT, maxValue_);
}

// ************************************************************************* //

