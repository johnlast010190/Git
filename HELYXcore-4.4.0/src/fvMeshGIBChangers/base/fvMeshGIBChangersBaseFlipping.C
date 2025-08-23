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

#include "base/fvMeshGIBChangersBase.H"
#include "meshes/meshTools/simpleVTKWriter.H"
#include "fvMesh/fvPatches/constraint/empty/emptyFvPatch.H"
#include "fvMesh/fvPatches/constraint/symmetry/symmetryFvPatch.H"
#include "fvMesh/fvPatches/constraint/wedge/wedgeFvPatch.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::fvMeshGIBChangersBase::flipCells
(
    boolList& interPoints,
    bool regionType
) const
{
    interPoints = false;
    boolList markInterFace(mesh().faces().size(), false);

    const labelList& faceZoneAddressing = fl();

    // Mark interface points
    forAll(faceZoneAddressing, i)
    {
        const face& f = mesh().faces()[faceZoneAddressing[i]];
        markInterFace[faceZoneAddressing[i]] = true;
        forAll(f, pointi)
        {
            interPoints[f[pointi]] = true;
        }
    }
    syncTools::syncFaceList(mesh(), markInterFace, orEqOp<bool>());
    syncTools::syncPointList(mesh(), interPoints, plusEqOp<bool>(), true);

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

    boolList cellId(mesh().cells().size(), false);
    boolList faceId(mesh().faces().size(), false);
    const labelList& reg = cRegion();
    forAll(reg, celli)
    {
        cellId[reg[celli]] = true;
    }
    forAll(faceId, facei)
    {
        if (facei < mesh().nInternalFaces())
        {
            faceId[facei] =
                !cellId[mesh().owner()[facei]]^cellId[mesh().neighbour()[facei]];
        }
        else
        {
            faceId[facei] = cellId[mesh().faceOwner()[facei]];
        }
    }
    syncTools::syncFaceList(mesh(), faceId, xxorEqOp());

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

    syncTools::syncFaceList(mesh(), markInterFace, orEqOp<bool>());

    singleCellSuddenPopCellRemoval(markInterFace);

    interPoints = false;
    DynamicList<label> dfl(mesh().faces().size());
    forAll(mesh().faces(), facei)
    {
        if (markInterFace[facei])
        {
            if (facei < mesh().nInternalFaces())
            {
                forAll(mesh().faces()[facei], facePointi)
                {
                    const label pointi = mesh().faces()[facei][facePointi];
                    interPoints[pointi] = true;
                }
                dfl.append(facei);
            }
            else
            {
                const label patchi = mesh().boundaryMesh().whichPatch(facei);
                const polyPatch& pp = mesh().boundary()[patchi].patch();
                if (pp.coupled() || includeWalls())
                {
                    forAll(mesh().faces()[facei], pfI)
                    {
                        const label pointi = mesh().faces()[facei][pfI];
                        interPoints[pointi] = true;
                    }
                    dfl.append(facei);
                }
            }
        }
    }

    syncTools::syncPointList(mesh(), interPoints, plusEqOp<bool>(), true);

    writeProblematicCells(interPoints);

    dfl.shrink();

    if (false)
    {
        const pointField& basePoints = *basePoints_;
        simpleVTKWriter
        (
            mesh().faces(),
            dfl,
            basePoints
        ).write
        (
            "fl0Flipped" + mesh().time().timeName() + ".vtk"
        );
    }

    delete flPtr_;
    flPtr_ = new labelList(dfl);

    deleteDemandDrivenData(cRegionPtr_);
    deleteDemandDrivenData(fmPtr_);
}


void Foam::fvMeshGIBChangersBase::singleCellSuddenPopCellRemoval
(
    boolList& markInterface
) const
{
    const boolList& markInterface0 = faceIndicator0();

    label nInternalFaces = mesh().nInternalFaces();
    forAll(mesh().cells(), celli)
    {
        const labelList& cFaces = mesh().cells()[celli];
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


Foam::boolList Foam::fvMeshGIBChangersBase::checkFlipCell
(
    const boolList& interPoints
) const
{
    const labelList& fll = fl();
    boolList flipCell(mesh().cells().size(), false);

    forAll(fll, fllI)
    {
        const label facei = fll[fllI];
        bool flip = true;
        if (facei < mesh().nInternalFaces())
        {
            const label own = mesh().owner()[facei];
            const label nei = mesh().neighbour()[facei];
            const labelList& cPointsOwn = mesh().cellPoints()[own];
            const labelList& cPointsNei = mesh().cellPoints()[nei];
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
            const label own = mesh().faceOwner()[facei];
            const labelList& cPointsOwn = mesh().cellPoints()[own];
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


Foam::labelList Foam::fvMeshGIBChangersBase::faceFlippingType
(
    const boolList& flipCell,
    const boolList& markInterFace
) const
{
    labelList faceFlipType(mesh().faces().size(), 0);
    forAll(flipCell, celli)
    {
        if (flipCell[celli])
        {
            const labelList& cFaces = mesh().cells()[celli];
            forAll(cFaces, cellFacei)
            {
                const label facei = cFaces[cellFacei];
                if (facei < mesh().nInternalFaces())
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
                        mesh().boundaryMesh().whichPatch(facei);
                    if
                    (
                        !isA<wedgeFvPatch>(mesh().boundary()[patchi])
                     && !isA<emptyFvPatch>(mesh().boundary()[patchi])
                     && !isA<symmetryFvPatch>(mesh().boundary()[patchi])
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

    syncTools::syncFaceList(mesh(), faceFlipType, plusEqOp<label>());

    return faceFlipType;
}


Foam::labelList Foam::fvMeshGIBChangersBase::cellFlippingNumber
(
    const boolList& flipCell,
    const boolList& markInterFace
) const
{
    labelList cellFlippingN(mesh().cells().size(), 0);
    forAll(flipCell, celli)
    {
        if (flipCell[celli])
        {
            const labelList& cFaces = mesh().cells()[celli];
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


void Foam::fvMeshGIBChangersBase::findFlippingFaces
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
            if (facei < mesh().nInternalFaces())
            {
                const label own = mesh().owner()[facei];
                const label nei = mesh().neighbour()[facei];
                const label cellFlipNon = cellFlipN[own];
                const label cellFlipNnb = cellFlipN[nei];

                if
                (
                    (mesh().cellPoints()[own].size() == 8)
                 && (mesh().cellPoints()[nei].size() == 8)
                )
                {
                    if ((cellFlipNnb == 3) && (cellFlipNon == 4))
                    {
                        const labelList& cFacesOwn = mesh().cells()[own];
                        forAll(cFacesOwn, cellFacei)
                        {
                            const label cellFacesI = cFacesOwn[cellFacei];
                            markInterFace[cellFacesI] =
                                !markInterFace[cellFacesI];
                            faceFlipType[cellFacesI] = 0;
                        }

                        const labelList& cFacesNei = mesh().cells()[nei];
                        forAll(cFacesNei, cellFacei)
                        {
                            faceFlipType[cFacesNei[cellFacei]] = 0;
                        }
                    }

                    if ((cellFlipNnb == 4) && (cellFlipNon == 3))
                    {
                        const labelList& cFacesNei = mesh().cells()[nei];
                        forAll(cFacesNei, cellFacei)
                        {
                            const label cellFacesI = cFacesNei[cellFacei];
                            markInterFace[cellFacesI] =
                                !markInterFace[cellFacesI];
                            faceFlipType[cellFacesI] = 0;
                        }

                        const labelList& cFacesOwn = mesh().cells()[own];
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


void Foam::fvMeshGIBChangersBase::resetFlipCellsWithFewInterfaces
(
    boolList& flipCell,
    const labelList& cellFlipN
) const
{
    forAll(flipCell, celli)
    {
        if (flipCell[celli])
        {
            if (mesh().cells()[celli].size() == 6)
            {
                flipCell[celli] = false;
            }
            else if (mesh().cells()[celli].size() == 5)
            {
                flipCell[celli] = false;
            }
        }
    }
}


void Foam::fvMeshGIBChangersBase::resetFlipHexCellsWith3Interfaces
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
                (mesh().cells()[celli].size() == 6)
             && (mesh().cellPoints()[celli].size() == 8)
            )
            {
                if (cellFlipN[celli] == 3)
                {
                    bool multipleflipCFould = false;
                    const labelList& cFaces = mesh().cells()[celli];
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
            if (mesh().cells()[celli].size() == 5)
            {
                if (cellFlipN[celli] == 2)
                {
                    bool multipleflipCFould = false;
                    const labelList& cFaces = mesh().cells()[celli];
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
