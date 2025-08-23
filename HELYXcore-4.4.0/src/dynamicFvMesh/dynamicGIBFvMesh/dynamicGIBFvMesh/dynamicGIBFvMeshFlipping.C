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
    (c) 2017 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "dynamicGIBFvMesh/dynamicGIBFvMesh/dynamicGIBFvMesh.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "fvMesh/fvPatches/constraint/wedge/wedgeFvPatch.H"
#include "fvMesh/fvPatches/constraint/empty/emptyFvPatch.H"
#include "fvMesh/fvPatches/constraint/symmetry/symmetryFvPatch.H"
#include "meshes/meshTools/simpleVTKWriter.H"


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::dynamicGIBFvMesh::flipCells
(
    boolList& interPoints,
    bool regionType
) const
{
    interPoints = false;
    boolList markInterFace(this->faces().size(), false);

    const labelList& faceZoneAddressing = fl();

    // Mark interface points
    forAll(faceZoneAddressing, i)
    {
        const face& f = this->faces()[faceZoneAddressing[i]];
        markInterFace[faceZoneAddressing[i]] = true;
        forAll(f, pointi)
        {
            interPoints[f[pointi]] = true;
        }
    }
    syncTools::syncFaceList(*this, markInterFace, orEqOp<bool>());
    syncTools::syncPointList(*this, interPoints, plusEqOp<bool>(), true);

    // Mark all the cells in which all the points are interface points
    boolList flipCell = checkFlipCell(interPoints);

/*
    forAll(flipCell, celli)
    {
        if (flipCell[celli])
        {
            if
            (
                (reg[celli] == 0 && regionType == false)
             || (reg[celli] == 1 && regionType == true)
            )
            {
               flipCell[cI] = false;
            }
        }
    }
*/

    labelList cellFlipN = cellFlippingNumber(flipCell, markInterFace);

    if (!pyrPrismFlip_)
    {
        resetFlipCellsWithFewInterfaces(flipCell, cellFlipN);
    }

    labelList faceFlipType = faceFlippingType(flipCell, markInterFace);

/*
    findFlippingFaces
    (
        flipCell,
        markInterFace,
        interPoints,
        cellFlipN,
        faceFlipType
    );

   resetFlipHexCellsWith3Interfaces(flipCell, cellFlipN, faceFlipType);
*/

    boolList cellId(this->cells().size(), false);
    boolList faceId(this->faces().size(), false);
    const labelList& reg = cRegion();
    forAll(reg, celli)
    {
        cellId[reg[celli]] = true;
    }
    forAll(faceId, facei)
    {
        if (facei < nInternalFaces())
        {
            faceId[facei] =
                !cellId[this->owner()[facei]]^cellId[this->neighbour()[facei]];
        }
        else
        {
            faceId[facei] = cellId[this->faceOwner()[facei]];
        }
    }
    syncTools::syncFaceList(*this, faceId, xxorEqOp());

    // Snap - unsnap cell based on faces
    forAll(markInterFace, facei)
    {
        if (faceFlipType[facei] == 1)
        {
            markInterFace[facei] = true;
        }
        else if (faceFlipType[facei] == -1)
        {
            markInterFace[facei] = false;
        }
        else if (faceFlipType[facei] == 2)
        {
            if (faceId[facei])
            {
                markInterFace[facei] = false;
            }
            else
            {
                markInterFace[facei] = true;
            }
        }
        else if (faceFlipType[facei] == -2)
        {
            markInterFace[facei] = true;
        }
    }

    syncTools::syncFaceList(*this, markInterFace, orEqOp<bool>());

    singleCellSuddenPopCellRemoval(markInterFace);

    interPoints = false;
    DynamicList<label> dfl(this->faces().size());
    forAll(this->faces(), facei)
    {
        if (markInterFace[facei])
        {
            if (facei < this->nInternalFaces())
            {
                forAll(this->faces()[facei], facePointi)
                {
                    const label pointi = this->faces()[facei][facePointi];
                    interPoints[pointi] = true;
                }
                dfl.append(facei);
            }
            else
            {
                const label patchi = this->boundaryMesh().whichPatch(facei);
                const polyPatch& pp = this->boundary()[patchi].patch();
                if (pp.coupled() || includeWalls())
                {
                    forAll(this->faces()[facei], pfI)
                    {
                        const label pointi = this->faces()[facei][pfI];
                        interPoints[pointi] = true;
                    }
                    dfl.append(facei);
                }
            }
        }
    }

    syncTools::syncPointList(*this, interPoints, plusEqOp<bool>(), true);

    writeProblematicCells(interPoints);

    dfl.shrink();


    if (false)
    {
        const pointField& basePoints = *basePoints_;
        simpleVTKWriter
        (
            this->faces(),
            dfl,
            basePoints
        ).write
        (
            "fl0Flipped" + this->time().timeName() + ".vtk"
        );
    }

    delete flPtr_;
    flPtr_ = new labelList(dfl);

    deleteDemandDrivenData(cRegionPtr_);
    deleteDemandDrivenData(fmPtr_);

}


void Foam::dynamicGIBFvMesh::singleCellSuddenPopCellRemoval
(
    boolList& markInterface
) const
{
    const boolList& markInterface0 = faceIndicator0();

    label nInternalFaces = this->nInternalFaces();
    forAll(this->cells(), celli)
    {
        const labelList& cFaces = this->cells()[celli];
        bool unsnap = true;
        forAll(cFaces, cellFacei)
        {
            const label facei = cFaces[cellFacei];
            if
            (
                markInterface0[facei]
            ||  !markInterface[facei]
            ||  !(facei < nInternalFaces)
            )
            {
                unsnap = false;
            }
        }
        if (unsnap)
        {
            forAll(cFaces, cellFacei)
            {
                markInterface[cFaces[cellFacei]] = false;
            }
        }
    }
}


Foam::boolList Foam::dynamicGIBFvMesh::checkFlipCell
(
    const boolList& interPoints
) const
{
    const labelList& fll = fl();
    boolList flipCell(this->cells().size(), false);

    forAll(fll, fllI)
    {
        const label facei = fll[fllI];
        bool flip = true;
        if (facei < this->nInternalFaces())
        {
            const label own = this->owner()[facei];
            const label nei = this->neighbour()[facei];
            const labelList& cPointsOwn = this->cellPoints()[own];
            const labelList& cPointsNei = this->cellPoints()[nei];
            forAll(cPointsOwn, cellPointi)
            {
                const label pointi = cPointsOwn[cellPointi];
                if (!interPoints[pointi])
                {
                    flip = false;
                }
            }
            flipCell[own] = flip;

            flip = true;
            forAll(cPointsNei, cellPointi)
            {
                const label pointi = cPointsNei[cellPointi];
                if (!interPoints[pointi])
                {
                   flip = false;
                }
            }
            flipCell[nei] = flip;
        }
        else
        {
            flip = true;
            const label own = this->faceOwner()[facei];
            const labelList& cPointsOwn = this->cellPoints()[own];
            forAll(cPointsOwn, cellPointi)
            {
                const label pointi = cPointsOwn[cellPointi];
                if (!interPoints[pointi])
                {
                    flip = false;
                }
            }
            flipCell[own] = flip;
        }
    }
    return flipCell;
}


Foam::labelList Foam::dynamicGIBFvMesh::faceFlippingType
(
    const boolList& flipCell,
    const boolList& markInterFace
) const
{
    labelList faceFlipType(this->faces().size(), 0);
    forAll(flipCell, celli)
    {
        if (flipCell[celli])
        {
            const labelList& cFaces = this->cells()[celli];
            forAll(cFaces, cellFacei)
            {
                const label facei = cFaces[cellFacei];
                if (facei < this->nInternalFaces())
                {
                    if (!markInterFace[facei])
                    {
                        faceFlipType[facei] += 1;
                    }
                    else
                    {
                        faceFlipType[facei] -= 1;
                    }
                }
                else
                {
                    const label patchi =
                        this->boundaryMesh().whichPatch(facei);
                    if
                    (
                        !isA<wedgeFvPatch>(this->boundary()[patchi])
                     && !isA<emptyFvPatch>(this->boundary()[patchi])
                     && !isA<symmetryFvPatch>(this->boundary()[patchi])
                    )
                    {
                        if (!markInterFace[facei])
                        {
                            faceFlipType[facei] += 1;
                        }
                        else
                        {
                            faceFlipType[facei] -= 1;
                        }
                    }
                }
            }
        }
    }

    syncTools::syncFaceList(*this, faceFlipType, plusEqOp<label>());
    return faceFlipType;
}


Foam::labelList Foam::dynamicGIBFvMesh::cellFlippingNumber
(
    const boolList& flipCell,
    const boolList& markInterFace
) const
{
    labelList cellFlippingN(this->cells().size(), 0);
    forAll(flipCell, celli)
    {
        if (flipCell[celli])
        {
            const labelList& cFaces = this->cells()[celli];
            forAll(cFaces, cellFacei)
            {
                const label facei = cFaces[cellFacei];
                if (markInterFace[facei])
                {
                    cellFlippingN[celli] += 1;
                }
            }
        }
    }
    return cellFlippingN;
}


void Foam::dynamicGIBFvMesh::findFlippingFaces
(
    const boolList& flipCell,
    boolList& markInterFace,
    const boolList& interPoints,
    const labelList& cellFlipN,
    labelList& faceFlipType
) const
{
    forAll(faceFlipType, facei)
    {
        if (faceFlipType[facei] == -2)
        {
            if (facei < this->nInternalFaces())
            {
                const label own = this->owner()[facei];
                const label nei = this->neighbour()[facei];
                const label cellFlipNon = cellFlipN[own];
                const label cellFlipNnb = cellFlipN[nei];

                if
                (
                    (this->cellPoints()[own].size() == 8)
                 && (this->cellPoints()[nei].size() == 8)
                )
                {
                    if ((cellFlipNnb == 3) && (cellFlipNon == 4))
                    {
                        const labelList& cFacesOwn = this->cells()[own];
                        forAll(cFacesOwn, cellFacei)
                        {
                            const label cellFacesI = cFacesOwn[cellFacei];
                            markInterFace[cellFacesI] =
                                !markInterFace[cellFacesI];
                            faceFlipType[cellFacesI] = 0;
                        }

                        const labelList& cFacesNei = this->cells()[nei];
                        forAll(cFacesNei, cellFacei)
                        {
                            faceFlipType[cFacesNei[cellFacei]] = 0;
                        }
                    }


                    if ((cellFlipNnb == 4) && (cellFlipNon == 3))
                    {
                        const labelList& cFacesNei = this->cells()[nei];
                        forAll(cFacesNei, cellFacei)
                        {
                            const label cellFacesI = cFacesNei[cellFacei];
                            markInterFace[cellFacesI] =
                                !markInterFace[cellFacesI];
                            faceFlipType[cellFacesI] = 0;
                        }

                        const labelList& cFacesOwn = this->cells()[own];
                        forAll(cFacesOwn, cellFacei)
                        {
                            faceFlipType[cFacesOwn[cellFacei]] = 0;
                        }
                    }
                }
            }
        }
    }
}


void Foam::dynamicGIBFvMesh::resetFlipCellsWithFewInterfaces
(
    boolList& flipCell,
    const labelList& cellFlipN
) const
{
    forAll(flipCell, celli)
    {
        if (flipCell[celli])
        {
            if (this->cells()[celli].size() == 6)
            {
                flipCell[celli] = false;
            }
            else if (this->cells()[celli].size() == 5)
            {
                flipCell[celli] = false;
            }
        }
    }
}


void Foam::dynamicGIBFvMesh::resetFlipHexCellsWith3Interfaces
(
    boolList& flipCell,
    const labelList& cellFlipN,
    labelList& faceFlipType
) const
{
    forAll(flipCell, celli)
    {
        if (flipCell[celli])
        {
            if
            (
                (this->cells()[celli].size() == 6)
             && (this->cellPoints()[celli].size() == 8)
            )
            {
                if (cellFlipN[celli] == 3)
                {
                    bool multipleflipCFould = false;
                    const labelList& cFaces = this->cells()[celli];
                    forAll(cFaces, cellFacei)
                    {
                        const label facei = cFaces[cellFacei];
                        if
                        (
                            (faceFlipType[facei] == 2)
                         || (faceFlipType[facei] == -2)
                        )
                        {
                            multipleflipCFould = true;
                        }
                    }
                    if (!multipleflipCFould)
                    {
                        flipCell[celli] = false;
                        forAll(cFaces, cellFacei)
                        {
                            faceFlipType[cFaces[cellFacei]] = 0;
                        }
                    }
                }
            }
            if (this->cells()[celli].size() == 5)
            {
                if (cellFlipN[celli] == 2)
                {
                    bool multipleflipCFould = false;
                    const labelList& cFaces = this->cells()[celli];
                    forAll(cFaces, cellFacei)
                    {
                        const label facei = cFaces[cellFacei];
                        if
                        (
                            (faceFlipType[facei] == 2)
                         || (faceFlipType[facei] == -2)
                        )
                        {
                            multipleflipCFould = true;
                        }
                    }
                    if (!multipleflipCFould)
                    {
                        flipCell[celli] = false;
                        forAll(cFaces, cellFacei)
                        {
                            faceFlipType[cFaces[cellFacei]] = 0;
                        }
                    }
                }
            }
        }
    }
}

// ************************************************************************* //
