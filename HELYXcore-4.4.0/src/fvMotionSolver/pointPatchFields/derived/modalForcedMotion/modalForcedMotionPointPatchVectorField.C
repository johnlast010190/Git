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

#include "pointPatchFields/derived/modalForcedMotion/modalForcedMotionPointPatchVectorField.H"
#include "fields/pointPatchFields/pointPatchField/pointPatchFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

modalForcedMotionPointPatchVectorField::modalForcedMotionPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    bufferZoneDisplacementBasePointPatchVectorField(p, iF),
    radialDisplacement_(),
    amplification_()
{}


modalForcedMotionPointPatchVectorField::modalForcedMotionPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    bufferZoneDisplacementBasePointPatchVectorField(p, iF, dict),
    origin_(dict.lookupOrDefault<vector>("origin", vector::zero)),
    axis_(dict.lookupOrDefault<vector>("axis", vector(0,0,1))),
    radialDisplacement_(Function1<vector>::New("radialDisplacement", dict)),
    amplification_(Function1<scalar>::New("amplification", dict))
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }
}


modalForcedMotionPointPatchVectorField::modalForcedMotionPointPatchVectorField
(
    const modalForcedMotionPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    bufferZoneDisplacementBasePointPatchVectorField(ptf, p, iF, mapper),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    radialDisplacement_(ptf.radialDisplacement_, false),
    amplification_(ptf.amplification_, false)
{}


modalForcedMotionPointPatchVectorField::modalForcedMotionPointPatchVectorField
(
    const modalForcedMotionPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    bufferZoneDisplacementBasePointPatchVectorField(ptf, iF),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    radialDisplacement_(ptf.radialDisplacement_, false),
    amplification_(ptf.amplification_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> modalForcedMotionPointPatchVectorField::getRadialCoordinates
(
    const vectorField& oldPoints
)
{
    tmp<scalarField> tradialCoordinates
    (
        new scalarField(oldPoints.size(), 0.0)
    );
    scalarField& radialCoordinates = tradialCoordinates.ref();

    // distance meansured in axis direction
    radialCoordinates = (oldPoints - origin_) & axis_;

    return tradialCoordinates;
}


void modalForcedMotionPointPatchVectorField::calcMaxDisplacement
(
    const vectorField& oldPoints,
    vectorField& maxDisplacement
)
{
    scalarField radialCoordinates(getRadialCoordinates(oldPoints));

    forAll(radialCoordinates, pI)
    {
        maxDisplacement[pI] =
            radialDisplacement_->value(radialCoordinates[pI]);
    }
}


tmp<vectorField> modalForcedMotionPointPatchVectorField::getDisplacement
(
    const vectorField& oldPoints
)
{
    const polyMesh& mesh = this->internalField().mesh()();
    const scalar t = mesh.time().value();

    vectorField maxDisplacement(oldPoints.size(), vector::zero);
    calcMaxDisplacement(oldPoints, maxDisplacement);

    return maxDisplacement * amplification_->value(t);
}


void modalForcedMotionPointPatchVectorField::write(Ostream& os) const
{
    bufferZoneDisplacementBasePointPatchVectorField::write(os);

    os.writeEntryIfDifferent<vector>
    (
        "origin",
        vector(0, 0, 0),
        origin_
    );
    os.writeEntryIfDifferent<vector>
    (
        "axis",
        vector(0, 0, 1),
        axis_
    );

    radialDisplacement_->writeData(os);
    amplification_->writeData(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    modalForcedMotionPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
