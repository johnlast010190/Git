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

#include "pointPatchFields/derived/bufferZoneDisplacementBase/bufferZoneDisplacementBasePointPatchVectorField.H"
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

bufferZoneDisplacementBasePointPatchVectorField::
bufferZoneDisplacementBasePointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    p0_(p.localPoints()),
    bufferSize_(0.0)
{}


bufferZoneDisplacementBasePointPatchVectorField::
bufferZoneDisplacementBasePointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    bufferSize_(dict.lookup<scalar>("bufferSize"))
{
    calculateBufferCellsAndPoints();

    //if (!dict.found("value"))
    //{
    //    updateCoeffs();
    //}

    if (dict.found("p0"))
    {
        p0_ = vectorField("p0", dict , p.size());
    }
    else
    {
        p0_ = p.localPoints();
    }
}


bufferZoneDisplacementBasePointPatchVectorField::
bufferZoneDisplacementBasePointPatchVectorField
(
    const bufferZoneDisplacementBasePointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    p0_(mapper(ptf.p0_)),
    bufferSize_(ptf.bufferSize_),
    bufferCells_(ptf.bufferCells_),
    bufferPoints_(ptf.bufferPoints_)
{}


bufferZoneDisplacementBasePointPatchVectorField::
bufferZoneDisplacementBasePointPatchVectorField
(
    const bufferZoneDisplacementBasePointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    p0_(ptf.p0_),
    bufferSize_(ptf.bufferSize_),
    bufferCells_(ptf.bufferCells_),
    bufferPoints_(ptf.bufferPoints_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void bufferZoneDisplacementBasePointPatchVectorField::autoMap
(
    const pointPatchFieldMapper& m
)
{
    fixedValuePointPatchField<vector>::autoMap(m);

    m(p0_, p0_);
}


void bufferZoneDisplacementBasePointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    const bufferZoneDisplacementBasePointPatchVectorField& aODptf =
        refCast<const bufferZoneDisplacementBasePointPatchVectorField>(ptf);

    fixedValuePointPatchField<vector>::rmap(aODptf, addr);

    p0_.rmap(aODptf.p0_, addr);
}


void bufferZoneDisplacementBasePointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    vectorField::operator=(this->getDisplacement(p0_));
    fixedValuePointPatchField<vector>::updateCoeffs();
}


void bufferZoneDisplacementBasePointPatchVectorField::manipulateMatrix
(
    const vectorField& oldPoints,
    labelList& bufferCells,
    vectorField& displacement
)
{
    const polyMesh& mesh = this->internalField().mesh()();
    bufferCells = bufferCells_;

    // Compute old cell centres given old point locations
    vectorField oldCellCentres(bufferCells.size(), vector::zero);
    forAll(bufferCells, i)
    {
        label cellI = bufferCells[i];
        cell c = mesh.cells()[cellI];
        oldCellCentres[i] = c.centre(oldPoints, mesh.faces());
    }

    // Compute cell centre displacement in buffer zone
    displacement = this->getDisplacement(oldCellCentres);
}


void bufferZoneDisplacementBasePointPatchVectorField::setField
(
    const vectorField& oldPoints,
    vectorField& displacement,
    bool setAll
)
{
    if (setAll)
    {
        displacement = this->getDisplacement(oldPoints);
    }
    else
    {
        // Create pointField in buffer zone
        pointField points(oldPoints, bufferPoints_);

        // Evaluate displacement only in buffer zone, then add entries
        UIndirectList<point>(displacement, bufferPoints_) =
            this->getDisplacement(points);
    }
}


void Foam::bufferZoneDisplacementBasePointPatchVectorField::calculateBufferCellsAndPoints()
{
    const polyMesh& mesh = this->internalField().mesh()();
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    label patchI = bMesh.findPatchID(this->patch().name());

    if (patchI >= 0 && bufferSize_ > SMALL)
    {
        labelHashSet patchSet(1);
        patchSet.insert(patchI);

        scalarField y(mesh.nCells());
        patchDistWave::wave<wallPoint>(mesh, patchSet, y, false);

        boolList isBufferPt(mesh.nPoints(), false);

        DynamicList<label> cells(mesh.nCells());
        DynamicList<label> points(mesh.nCells());
        forAll(y, cellI)
        {
            if (y[cellI] < bufferSize_)
            {
                cells.append(cellI);
                cell c = mesh.cells()[cellI];
                labelList meshPoints = c.labels(mesh.faces());
                forAll(meshPoints, ptI)
                {
                    isBufferPt[meshPoints[ptI]] = true;
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

tmp<vectorField> bufferZoneDisplacementBasePointPatchVectorField::getDisplacement
(
    const vectorField& oldPoints
)
{
    // NotImplemented;
    return tmp<vectorField>
    (
        new vectorField(oldPoints.size(), vector::zero)
    );
}


void bufferZoneDisplacementBasePointPatchVectorField::write
(
    Ostream& os
) const
{
    pointPatchField<vector>::write(os);

    os.writeEntry("bufferSize", bufferSize_);
    p0_.writeEntry("p0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    bufferZoneDisplacementBasePointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
