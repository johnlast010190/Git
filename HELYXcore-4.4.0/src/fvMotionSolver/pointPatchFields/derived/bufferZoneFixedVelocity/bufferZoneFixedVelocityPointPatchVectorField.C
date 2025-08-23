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

#include "pointPatchFields/derived/bufferZoneFixedVelocity/bufferZoneFixedVelocityPointPatchVectorField.H"
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bufferZoneFixedVelocityPointPatchVectorField::
bufferZoneFixedVelocityPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    radialVelocity_(),
    linearVelocity_(),
    CofG_(vector::zero),
    bufferSize_(0.0)
{}


Foam::bufferZoneFixedVelocityPointPatchVectorField::
bufferZoneFixedVelocityPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    radialVelocity_(Function1<vector>::New("radialVelocity", dict)),
    linearVelocity_(Function1<vector>::New("linearVelocity", dict)),
    CofG_(dict.lookup("CofG")),
    bufferSize_(dict.lookup<scalar>("bufferSize"))
{
    calculateBufferCellsAndPoints();

    if (!dict.found("value"))
    {
        bufferZoneFixedVelocityPointPatchVectorField::updateCoeffs();
    }
}


Foam::bufferZoneFixedVelocityPointPatchVectorField::
bufferZoneFixedVelocityPointPatchVectorField
(
    const bufferZoneFixedVelocityPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    radialVelocity_(ptf.radialVelocity_().clone().ptr()),
    linearVelocity_(ptf.linearVelocity_().clone().ptr()),
    CofG_(ptf.CofG_),
    bufferSize_(ptf.bufferSize_)
{}


Foam::bufferZoneFixedVelocityPointPatchVectorField::
bufferZoneFixedVelocityPointPatchVectorField
(
    const bufferZoneFixedVelocityPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    radialVelocity_(ptf.radialVelocity_().clone().ptr()),
    linearVelocity_(ptf.linearVelocity_().clone().ptr()),
    CofG_(ptf.CofG_),
    bufferSize_(ptf.bufferSize_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::bufferZoneFixedVelocityPointPatchVectorField::autoMap
(
    const pointPatchFieldMapper& m
)
{
    fixedValuePointPatchField<vector>::autoMap(m);
}


void Foam::bufferZoneFixedVelocityPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    const bufferZoneFixedVelocityPointPatchVectorField& aODptf =
        refCast<const bufferZoneFixedVelocityPointPatchVectorField>(ptf);

    fixedValuePointPatchField<vector>::rmap(aODptf, addr);
}


void Foam::bufferZoneFixedVelocityPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();

    vector transVec = linearVelocity_->value(t.timeOutputValue());
    vector eulerAngles = radialVelocity_->value(t.timeOutputValue());

    eulerAngles *= Foam::constant::mathematical::pi*t.deltaTValue()/180.0;

    quaternion R(quaternion::XYZ, eulerAngles);
    septernion TR(septernion(-CofG_ + -transVec)*R*septernion(CofG_));

    vectorField lp = this->patch().localPoints();
    this->forceAssign((transformPoints(TR, lp) - lp)/t.deltaTValue());
    fixedValuePointPatchField<vector>::updateCoeffs();
}


void Foam::bufferZoneFixedVelocityPointPatchVectorField::manipulateMatrix
(
    labelList& bufferCells,
    vectorField& velocity
)
{
    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();

    vector transVec = linearVelocity_->value(t.timeOutputValue());
    vector eulerAngles = radialVelocity_->value(t.timeOutputValue());

    eulerAngles *= Foam::constant::mathematical::pi*t.deltaTValue()/180.0;

    quaternion R(quaternion::XYZ, eulerAngles);
    septernion TR(septernion(-CofG_ - transVec)*R*septernion(CofG_));

    bufferCells = bufferCells_;

    vectorField bufferCellCentres(bufferCells.size(), vector::zero);
    forAll(bufferCells, i)
    {
        label cellI = bufferCells[i];
        bufferCellCentres[i] = mesh.cellCentres()[cellI];
    }

    velocity =
        (
            transformPoints(TR, bufferCellCentres)
          - bufferCellCentres
        )/t.deltaTValue();
}


void Foam::bufferZoneFixedVelocityPointPatchVectorField::setField
(
    Field<vector>& displacement
)
{
    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();

    vector transVec = linearVelocity_->value(t.timeOutputValue());
    vector eulerAngles = radialVelocity_->value(t.timeOutputValue());

    eulerAngles *= Foam::constant::mathematical::pi*t.deltaTValue()/180.0;

    quaternion R(quaternion::XYZ, eulerAngles);
    septernion TR(septernion(-CofG_ + -transVec)*R*septernion(CofG_));

    UIndirectList<point>(displacement, bufferPoints_) =
        (
            transformPoints(TR, pointField(mesh.points(), bufferPoints_))
          - pointField(mesh.points(), bufferPoints_)
        )/t.deltaTValue();
}


void Foam::bufferZoneFixedVelocityPointPatchVectorField::
calculateBufferCellsAndPoints()
{
    const polyMesh& mesh = this->internalField().mesh()();
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    label patchI = bMesh.findPatchID(this->patch().name());

    if (patchI)
    {
        labelHashSet patchSet(1);
        patchSet.insert(patchI);

        scalarField y(mesh.nCells());
        patchDistWave::wave<wallPoint>(mesh, patchSet, y, false);

        boolList isBufferPt(mesh.nPoints(), true);

        DynamicList<label> cells(mesh.nCells());
        DynamicList<label> points(mesh.nCells());
        forAll(y, cellI)
        {
            if (y[cellI] < bufferSize_)
            {
                cells.append(cellI);
            }
            else
            {
                cell c = mesh.cells()[cellI];
                labelList meshPoints = c.labels(mesh.faces());
                forAll(meshPoints, ptI)
                {
                    isBufferPt[meshPoints[ptI]] = false;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            isBufferPt,
            orEqOp<bool>(),
            false
         );

        forAll(isBufferPt, ptI)
        {
            if (isBufferPt[ptI])
            {
               points.append(ptI);
            }
        }
        cells.shrink();
        points.shrink();

        bufferCells_ = cells;
        bufferPoints_ = points;
    }
    else
    {
        bufferCells_ = labelList();
        bufferPoints_ = labelList();
    }
}


void Foam::bufferZoneFixedVelocityPointPatchVectorField::write
(
    Ostream& os
) const
{
    pointPatchField<vector>::write(os);
    radialVelocity_->writeData(os);
    linearVelocity_->writeData(os);
    os.writeEntry("CofG", CofG_);
    os.writeEntry("bufferSize", bufferSize_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
    makePointPatchTypeField
    (
        pointPatchVectorField,
        bufferZoneFixedVelocityPointPatchVectorField
    );
}

// ************************************************************************* //
