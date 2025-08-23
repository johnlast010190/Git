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

#include "pointPatchFields/derived/sedimentation/sedimentationPointPatchVectorField.H"
#include "fields/pointPatchFields/pointPatchField/pointPatchFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "cfdTools/general/include/fvCFD.H"
#include "interpolations/primitivePatchInterpolation/primitivePatchInterpolation.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sedimentationPointPatchVectorField::sedimentationPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    constrainDisplacement_(false),
    sedimentationSurfaceConstrainName_(word::null),
    sedimentationSurfaceConstrain_(nullptr),
    normal_(0, 0, 1),
    isDisplacement_(false)
{}


Foam::sedimentationPointPatchVectorField::sedimentationPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
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
    isDisplacement_(dict.lookupOrDefault<Switch>("isDisplacement", false))
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
            << sedimentationPointPatchVectorField::typeName
            << " didn't find gravity vector g or upDirection"
            << " in the dictionary. Using default value (0, 0, 1)."
            << nl;
    }
    normal_ = normalised(normal_);
}


Foam::sedimentationPointPatchVectorField::sedimentationPointPatchVectorField
(
    const sedimentationPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    constrainDisplacement_(ptf.constrainDisplacement_),
    sedimentationSurfaceConstrainName_(ptf.sedimentationSurfaceConstrainName_),
    sedimentationSurfaceConstrain_(ptf.sedimentationSurfaceConstrain_),
    normal_(ptf.normal_),
    isDisplacement_(ptf.isDisplacement_)
{}


Foam::sedimentationPointPatchVectorField::sedimentationPointPatchVectorField
(
    const sedimentationPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    constrainDisplacement_(ptf.constrainDisplacement_),
    sedimentationSurfaceConstrainName_(ptf.sedimentationSurfaceConstrainName_),
    sedimentationSurfaceConstrain_(ptf.sedimentationSurfaceConstrain_),
    normal_(ptf.normal_),
    isDisplacement_(ptf.isDisplacement_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sedimentationPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const objectRegistry& obr = this->db();
    const Time& runTime = obr.time();
    const scalar deltaT = runTime.deltaTValue();

    const label patchi = this->patch().index();
    const volScalarField& deltaH = obr.lookupObject<volScalarField>("deltaH");
    const scalarField& deltaHp = deltaH.boundaryField()[patchi];

    // Interpolation from face center to node
    const polyPatch& p = deltaH.mesh().boundaryMesh()[patchi];
    const pointField& points = p.localPoints();
    primitivePatch pPatch(SubList<face>(p.localFaces(), p.size()), points);
    tmp<scalarField> deltaHpoints =
        primitivePatchInterpolation(pPatch).faceToPointInterpolate(deltaHp);

    // Initialise the velocity to zero
    tmp<vectorField> tPointU(normal_*deltaHpoints()/deltaT);
    vectorField& pointU = tPointU.ref();

    // Limiting sedimentation by the surface beyond which the mesh can't move
    // anymore
    if (constrainDisplacement_)
    {
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
            sedimentationSurfaceConstrain_().findRayNearest(points, -normal_)
        );

        forAll(*this, i)
        {
            const scalar deltaHMax = ((nearestPoints[i] - points[i])&normal_);
            const scalar deltaHcomputed = (pointU[i]*deltaT&normal_);

            // Limits only erosion
            if
            (
                nearestPoints[i] != vector(GREAT, GREAT, GREAT) // complete miss (intersection not found)
             && pointU[i] != vector(0, 0, 0) // Velocity is zero anyways no limit to velocity needed
             && deltaHcomputed < 0
             && deltaHcomputed < deltaHMax
            )
            {

                pointU[i] = (nearestPoints[i] - points[i])/(isDisplacement_ ? deltaT : 1.0);
            }
        }
    }

    this->forceAssign(tPointU);

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void Foam::sedimentationPointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);

    os.writeEntry("constrainDisplacement", constrainDisplacement_);
    os.writeEntryIfDifferent<word>
    (
        "surfaceConstrainName",
        word::null,
        sedimentationSurfaceConstrainName_
    );
    os.writeEntryIfDifferent<vector>("upDirection", vector(0, 0, 1), normal_);
    os.writeEntryIfDifferent<Switch>
    (
        "isDisplacement",
        false,
        isDisplacement_
    );

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePointPatchTypeField
    (
        pointPatchVectorField,
        sedimentationPointPatchVectorField
    );
}

// ************************************************************************* //
