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
    (c) 2011 OpenFOAM Foundation
    (c) 2013-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "pointPatchFields/derived/bufferZoneFixedDisplacement/bufferZoneFixedDisplacementPointPatchVectorField.H"
#include "fields/pointPatchFields/pointPatchField/pointPatchFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "patchDist/patchDistWave/patchDistWave.H"
#include "patchDist/wallPoint/wallPoint.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "meshes/polyMesh/syncTools/syncTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

bufferZoneFixedDisplacementPointPatchVectorField::
bufferZoneFixedDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    bufferZoneDisplacementBasePointPatchVectorField(p, iF),
    radialVelocity_(vector::zero),
    linearVelocity_(vector::zero),
    CofG_(vector::zero)
{}


bufferZoneFixedDisplacementPointPatchVectorField::
bufferZoneFixedDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    bufferZoneDisplacementBasePointPatchVectorField(p, iF, dict),
    radialVelocity_(dict.lookup("radialVelocity")),
    linearVelocity_(dict.lookup("linearVelocity")),
    CofG_(dict.lookup("CofG"))
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }
}


bufferZoneFixedDisplacementPointPatchVectorField::
bufferZoneFixedDisplacementPointPatchVectorField
(
    const bufferZoneFixedDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    bufferZoneDisplacementBasePointPatchVectorField(ptf, p, iF, mapper),
    radialVelocity_(ptf.radialVelocity_),
    linearVelocity_(ptf.linearVelocity_),
    CofG_(ptf.CofG_)
{}


bufferZoneFixedDisplacementPointPatchVectorField::
bufferZoneFixedDisplacementPointPatchVectorField
(
    const bufferZoneFixedDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    bufferZoneDisplacementBasePointPatchVectorField(ptf, iF),
    radialVelocity_(ptf.radialVelocity_),
    linearVelocity_(ptf.linearVelocity_),
    CofG_(ptf.CofG_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> bufferZoneFixedDisplacementPointPatchVectorField::getDisplacement
(
    const vectorField& oldPoints
)
{
    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();

    vector transVec = linearVelocity_*t.value();

    vector eulerAngles = radialVelocity_*t.value();
    eulerAngles *= Foam::constant::mathematical::pi/180.0;

    quaternion R(quaternion::XYZ, eulerAngles);
    septernion TR(septernion(-CofG_ + -transVec)*R*septernion(CofG_));

    return tmp<vectorField>
    (
        new vectorField(transformPoints(TR, oldPoints) - oldPoints)
    );
}


void bufferZoneFixedDisplacementPointPatchVectorField::write
(
    Ostream& os
) const
{
    bufferZoneDisplacementBasePointPatchVectorField::write(os);
    os.writeEntry("radialVelocity", radialVelocity_);
    os.writeEntry("linearVelocity", linearVelocity_);
    os.writeEntry("CofG", CofG_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    bufferZoneFixedDisplacementPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
