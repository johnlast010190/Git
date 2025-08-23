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
    (c) 2013-2018 OpenFOAM Foundation
    (c) 2020 OpenCFD Ltd.
    (c) 2021 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "AMIInterpolation/AMIInterpolation/AMIMethod/partialFaceAreaWeightAMI/partialFaceAreaWeightAMI.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "AMIInterpolation/AMIInterpolation/AMIInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(partialFaceAreaWeightAMI, 0);
    addToRunTimeSelectionTable(AMIMethod, partialFaceAreaWeightAMI, components);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::partialFaceAreaWeightAMI::setNextFaces
(
    label& startSeedI,
    label& srcFacei,
    label& tgtFacei,
    const boolList& mapFlag,
    labelList& seedFaces,
    const DynamicList<label>& visitedFaces,
    const bool errorOnNotFound
) const
{
    faceAreaWeightAMI::setNextFaces
    (
        startSeedI,
        srcFacei,
        tgtFacei,
        mapFlag,
        seedFaces,
        visitedFaces,
        false // no error on not found
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::partialFaceAreaWeightAMI::
partialFaceAreaWeightAMI
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch,
    const scalarField& srcMagSf,
    const scalarField& tgtMagSf,
    const faceAreaIntersect::triangulationMode& triMode,
    const scalar& cosMatchAngle,
    const bool reverseTarget,
    const bool requireMatch
)
:
    faceAreaWeightAMI
    (
        srcPatch,
        tgtPatch,
        srcMagSf,
        tgtMagSf,
        triMode,
        cosMatchAngle,
        reverseTarget,
        requireMatch
    )
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::partialFaceAreaWeightAMI::~partialFaceAreaWeightAMI()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::partialFaceAreaWeightAMI::conformal() const
{
    return false;
}


void Foam::partialFaceAreaWeightAMI::calculate
(
    labelListList& srcAddress,
    scalarListList& srcWeights,
    scalarList& srcSymWeights,
    pointListList& srcCentroids,
    labelListList& tgtAddress,
    scalarListList& tgtWeights,
    scalarList& tgtSymWeights,
    label srcFacei,
    label tgtFacei
)
{
    bool ok =
        this->initialise
        (
            srcAddress,
            srcWeights,
            srcSymWeights,
            tgtAddress,
            tgtWeights,
            tgtSymWeights,
            srcFacei,
            tgtFacei
        );

    if (!ok)
    {
        return;
    }

    // temporary storage for addressing and weights
    List<DynamicList<label>> srcAddr(this->srcPatch_.size());
    List<DynamicList<scalar>> srcWght(srcAddr.size());
    List<scalar> srcSymWght(srcAddr.size(), Zero);
    List<DynamicList<label>> tgtAddr(this->tgtPatch_.size());
    List<DynamicList<scalar>> tgtWght(tgtAddr.size());
    List<scalar> tgtSymWght(tgtAddr.size(), Zero);
    List<DynamicList<point>> srcCtr(srcAddr.size());

    faceAreaWeightAMI::calcAddressing
    (
        srcAddr,
        srcWght,
        srcSymWght,
        srcCtr,
        tgtAddr,
        tgtWght,
        tgtSymWght,
        srcFacei,
        tgtFacei
    );

    // Check for badly covered faces
    if (this->restartUncoveredFaces())
    {
        this->restartUncoveredSourceFace
        (
            srcAddr,
            srcWght,
            srcSymWght,
            srcCtr,
            tgtAddr,
            tgtWght,
            tgtSymWght
        );
    }


    // transfer data to persistent storage
    forAll(srcAddr, i)
    {
        srcAddress[i].transfer(srcAddr[i]);
        srcWeights[i].transfer(srcWght[i]);
        srcSymWeights[i] = srcSymWght[i];
    }
    forAll(tgtAddr, i)
    {
        tgtAddress[i].transfer(tgtAddr[i]);
        tgtWeights[i].transfer(tgtWght[i]);
        tgtSymWeights[i] = tgtSymWght[i];
    }
}


// ************************************************************************* //
