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
    (c) 2011-2016 OpenFOAM Foundation

Description
    Efficient cell-centre calculation using face-addressing, face-centres and
    face-areas.

\*---------------------------------------------------------------------------*/

#include "meshes/primitiveMesh/primitiveMesh.H"
#include "primitives/bools/Switch/Switch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::primitiveMesh::calcCellCentresAndVols() const
{
    if (debug)
    {
        Pout<< "primitiveMesh::calcCellCentresAndVols() : "
            << "Calculating cell centres and cell volumes"
            << endl;
    }

    // It is an error to attempt to recalculate cellCentres
    // if the pointer is already set
    if (cellCentresPtr_ || cellVolumesPtr_)
    {
        FatalErrorInFunction
            << "Cell centres or cell volumes already calculated"
            << abort(FatalError);
    }

    // set the accumulated cell centre to zero vector
    cellCentresPtr_ = new vectorField(nCells());
    vectorField& cellCtrs = *cellCentresPtr_;

    // Initialise cell volumes to 0
    cellVolumesPtr_ = new scalarField(nCells());
    scalarField& cellVols = *cellVolumesPtr_;

    // Make centres and volumes
    makeCellCentresAndVols(faceCentres(), faceAreas(), cellCtrs, cellVols);

    if (debug)
    {
        Pout<< "primitiveMesh::calcCellCentresAndVols() : "
            << "Finished calculating cell centres and cell volumes"
            << endl;
    }
}


void Foam::primitiveMesh::cellCentreAndVol
(
    const label& celli,
    const cell& c,
    const labelList& owner,
    const vectorField& fCtrs,
    const vectorField& fAreas,
    vector& cellCentre,
    scalar& cellVol
)
{
    cellCentre = Zero;
    cellVol = 0.0;
#if defined(HELYX_SP)
    Field<Vector<double>> cellFaceCtrs(c.size());
    Field<Vector<double>> cellFaceAreas(c.size());

    forAll(c, fI)
    {
       label facei = c[fI];
       cellFaceCtrs[fI] = fCtrs[facei];
       if (owner[facei] != celli)
       {
           cellFaceAreas[fI] = -fAreas[facei];
       }
       else
       {
           cellFaceAreas[fI] = fAreas[facei];
       }
    }

    Vector<double> cc = 0.0;
    double cv = Zero;
    Vector<double> cEst = Zero;
    forAll(c, fI)
    {
        cEst += cellFaceCtrs[fI];
    }
    cEst /= c.size();
    double cellVolMag(0.0);

    forAll(c, fI)
    {
        Vector<double> relFC = cellFaceCtrs[fI] - cEst;

        // Calculate 3*face-pyramid volume
        double pyr3Vol = (cellFaceAreas[fI] & relFC);
        double pyr3VolMag = mag(pyr3Vol);

        // Calculate face-pyramid centre
        const Vector<double> pc = 0.75*relFC;

        // Accumulate volume-weighted face-pyramid centre
        cc += pyr3VolMag*pc;

        // Accumulate face-pyramid volume
        cv += pyr3Vol;
        cellVolMag += pyr3VolMag;
    }

    if (mag(cellVolMag) > VSMALL)
    {
       cc /= cellVolMag;
    }
    else
    {
       cc = cEst;
    }
    // relative correction to cc
    cc += cEst;
    cv /= 3.0;
    cellCentre = cc;
    cellVol = cv;
#else
    vector cEst = Zero;
    forAll(c, fI)
        cEst += fCtrs[c[fI]];
    cEst /= c.size();
    scalar cellVolMag(0.0);

    forAll(c, fI)
    {
        label facei = c[fI];
        point relFC = fCtrs[facei] - cEst;

        // Calculate 3*face-pyramid volume
        scalar pyr3Vol = (fAreas[facei] & relFC);

        if (owner[facei] != celli)
            pyr3Vol *= -1.0;

        scalar pyr3VolMag = mag(pyr3Vol);

        // Calculate face-pyramid centre
        const vector pc = 0.75*relFC;

        // Accumulate volume-weighted face-pyramid centre
        cellCentre += pyr3VolMag*pc;

        // Accumulate face-pyramid volume
        cellVol += pyr3Vol;
        cellVolMag += pyr3VolMag;
    }

    if (mag(cellVolMag) > VSMALL)
    {
       cellCentre /= cellVolMag;
    }
    else
    {
       cellCentre = cEst;
    }
    // relative correction to cc
    cellCentre += cEst;
    cellVol /= 3.0;
#endif
}


void Foam::primitiveMesh::makeCellCentresAndVols
(
    const vectorField& fCtrs,
    const vectorField& fAreas,
    vectorField& cellCtrs,
    scalarField& cellVols
) const
{
    const labelList& owner = faceOwner();
    const cellList& cs = cells();

    forAll(cs, celli)
    {
        const cell& c = cs[celli];
        cellCentreAndVol
        (
            celli,
            c,
            owner,
            fCtrs,
            fAreas,
            cellCtrs[celli],
            cellVols[celli]
        );
    }
}

// ************************************************************************* //
