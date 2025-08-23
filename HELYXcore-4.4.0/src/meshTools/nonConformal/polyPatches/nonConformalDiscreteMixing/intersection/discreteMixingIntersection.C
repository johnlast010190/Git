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
    (c) 2024-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "nonConformal/polyPatches/nonConformalDiscreteMixing/intersection/discreteMixingIntersection.H"
#include "nonConformal/polyPatches/nonConformalDiscreteMixing/nonConformalDiscreteMixingPolyPatch.H"
#include "fields/Fields/Field/SubField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(discreteMixingIntersection, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::discreteMixingIntersection::setAddressing()
{
    forAll(patchToPatchList_, intersecti)
    {
        // Source addressing
        const auto& srcTgtProcFaces =
            patchToPatchList_[intersecti].srcTgtProcFaces();
        forAll(srcTgtProcFaces, srci)
        {
            forAll(srcTgtProcFaces[srci], i)
            {
                srcAddress_[srci].append(srcTgtProcFaces[srci][i]);
            }
        }

        // Target addressing
        const auto& tgtSrcProcFaces =
            patchToPatchList_[intersecti].tgtSrcProcFaces();
        forAll(tgtSrcProcFaces, tgti)
        {
            forAll(tgtSrcProcFaces[tgti], i)
            {
                tgtAddress_[tgti].append(tgtSrcProcFaces[tgti][i]);
            }
        }
    }
}


void Foam::discreteMixingIntersection::calcWeights(const bool isSrc)
{
    const auto& addressing = isSrc ? srcAddress_ : tgtAddress_;
    auto& weights = isSrc ? srcWeights_ : tgtWeights_;

    List<List<DynamicList<scalar>>> ictWeights(patchToPatchList_.size());
    scalarList areaCoverage(addressing.size(), Zero);

    // For each original face, calculate the coupled weights and the total
    // coupling area across all rotational intersections.
    forAll(patchToPatchList_, intersecti)
    {
        ictWeights[intersecti].resize(addressing.size());

        const auto& couples =
            isSrc
          ? patchToPatchList_[intersecti].srcCouples()
          : patchToPatchList_[intersecti].tgtCouples();

        forAll(couples, facei)
        {
            ictWeights[intersecti][facei].resize(couples[facei].size());

            forAll(couples[facei], i)
            {
                const scalar a = mag(couples[facei][i].area);

                ictWeights[intersecti][facei][i] = a;
                areaCoverage[facei] += a;
            }
        }
    }

    // Normalise the weights
    forAll(patchToPatchList_, intersecti)
    {
        const auto& couples =
            isSrc
          ? patchToPatchList_[intersecti].srcCouples()
          : patchToPatchList_[intersecti].tgtCouples();

        forAll(couples, facei)
        {
            forAll(couples[facei], i)
            {
                ictWeights[intersecti][facei][i] /= areaCoverage[facei];
            }
        }
    }

    // Combine per-patch-face weights
    forAll(ictWeights, intersecti)
    {
        forAll(ictWeights[intersecti], facei)
        {
            forAll(ictWeights[intersecti][facei], i)
            {
                weights[facei].append(ictWeights[intersecti][facei][i]);
            }
        }
    }

    // Normalise the faces coverage and check for zero-weight faces
    const auto& areas = isSrc ? srcAreas_ : tgtAreas_;
    scalarList coverage = areaCoverage/areas;

    labelHashSet& zeroWeightFaces =
        isSrc ? srcZeroWeightFaces_ : tgtZeroWeightFaces_;
    forAll(coverage, facei)
    {
        if (coverage[facei] < SMALL)
        {
            zeroWeightFaces.insert(facei);
        }
    }

    const label nZeroWeight =
        returnReduce(zeroWeightFaces.size(), sumOp<label>());
    if (nZeroWeight)
    {
        Info<< nl;
        WarningInFunction
            << nZeroWeight << " faces detected with a zero coupling "
            << "weight for this non-conformal discrete mixing "
            << (isSrc ? "source" : "target") << " patch." << endl;
    }
}


void Foam::discreteMixingIntersection::setIntersectionMaps()
{
    forAll(patchToPatchList_, intersecti)
    {
        // Source map
        forAll(srcAddress_, facei)
        {
            srcIntersectionMap_[facei].append
            (
                labelList
                (
                    patchToPatchList_[intersecti].srcCouples()[facei].size(),
                    intersecti
                )
            );
        }

        // Target map
        forAll(tgtAddress_, facei)
        {
            tgtIntersectionMap_[facei].append
            (
                labelList
                (
                    patchToPatchList_[intersecti].tgtCouples()[facei].size(),
                    intersecti
                )
            );
        }
    }
}


void Foam::discreteMixingIntersection::agglomerate
(
    const scalarField& fineAreas,
    const List<List<remote>>& fineAddress,
    const List<scalarList>& fineWeights,
    const labelList& restrictAddress,
    const labelList& otherRestrictAddress,
    const bool isSrc
)
{
    const label coarseSize =
    (
        restrictAddress.size() ? max(restrictAddress) + 1 : 0
    );

    // Alias the patch geometrical quantities
    scalarField& areas = isSrc ? srcAreas_ : tgtAreas_;
    List<List<remote>>& address = isSrc ? srcAddress_ : tgtAddress_;
    List<scalarList>& weights = isSrc ? srcWeights_ : tgtWeights_;
    labelHashSet& zeroWeightFaces =
        isSrc ? srcZeroWeightFaces_ : tgtZeroWeightFaces_;

    // Agglomerate face areas
    areas.setSize(restrictAddress.size(), Zero);
    forAll(restrictAddress, facei)
    {
        const label coarseFacei = restrictAddress[facei];
        areas[coarseFacei] += fineAreas[facei];
    }

    // Agglomerate addressing and weights

    address.setSize(coarseSize);
    weights.setSize(coarseSize);

    // Gather restrict address on all processors
    List<labelList> procRestrictAddress(Pstream::nProcs());
    procRestrictAddress[Pstream::myProcNo()] = restrictAddress;
    Pstream::allGatherList(procRestrictAddress);

    // Also gather restrict address for the coupled side
    List<labelList> procOtherRestrictAddress(Pstream::nProcs());
    procOtherRestrictAddress[Pstream::myProcNo()] = otherRestrictAddress;
    Pstream::allGatherList(procOtherRestrictAddress);

    // Find remote index in a list
    auto findRemote = []
    (
        const List<remote>& newProcFaces,
        const label& coupledProci,
        const label& coarseCoupledFacei
    )
    {
        for (label i = 0; i < newProcFaces.size(); i++)
        {
            if
            (
                newProcFaces[i].proci == coupledProci
             && newProcFaces[i].elementi == coarseCoupledFacei
            )
            {
                return i;
            }
        }

        return label(-1);
    };

    forAll(fineAddress, facei)
    {
        const List<remote>& procFaces = fineAddress[facei];
        const scalarList& coupledWeights = fineWeights[facei];
        const scalar fineArea = fineAreas[facei];

        const label coarseFacei = restrictAddress[facei];

        List<remote>& newProcFaces = address[coarseFacei];
        scalarList& newWeights = weights[coarseFacei];

        forAll(procFaces, i)
        {
            const label coupledFacei = procFaces[i].elementi;
            const label coupledProci = procFaces[i].proci;

            const label coarseCoupledFacei =
                procOtherRestrictAddress[coupledProci][coupledFacei];

            const label index =
                findRemote(newProcFaces, coupledProci, coarseCoupledFacei);

            if (index == -1)
            {
                newProcFaces.append(remote(coupledProci, coarseCoupledFacei));
                newWeights.append(fineArea*coupledWeights[i]);
            }
            else
            {
                newWeights[index] += fineArea*coupledWeights[i];
            }
        }
    }

    // Normalise the weights
    forAll(weights, facei)
    {
        scalarList& w = weights[facei];
        const scalar totalWeight = sum(w);

        forAll(w, i)
        {
            w[i] /= totalWeight;
        }
    }

    // Reset zero-weight faces
    zeroWeightFaces.clear();
    scalarList coverage(coarseSize, Zero);
    forAll(weights, facei)
    {
        forAll(weights[facei], i)
        {
            coverage[facei] += weights[facei][i];
        }

        if (coverage[facei] < SMALL)
        {
            zeroWeightFaces.insert(facei);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::discreteMixingIntersection::discreteMixingIntersection
(
    const nonConformalDiscreteMixingPolyPatch& ownerPatch,
    const nonConformalDiscreteMixingPolyPatch& neighbourPatch
)
:
    patchToPatchList_(),
    transforms_(),
    srcAreas_(ownerPatch.origPatch().magFaceAreas()),
    tgtAreas_(neighbourPatch.origPatch().magFaceAreas()),
    srcAddress_(ownerPatch.origPatch().size()),
    tgtAddress_(neighbourPatch.origPatch().size()),
    srcWeights_(ownerPatch.origPatch().size()),
    tgtWeights_(neighbourPatch.origPatch().size()),
    srcZeroWeightFaces_(),
    tgtZeroWeightFaces_(),
    srcIntersectionMap_(ownerPatch.origPatch().size()),
    tgtIntersectionMap_(neighbourPatch.origPatch().size())
{
    // Define lists size as the total number of rotational intersections
    const label nIntersections =
        ownerPatch.patchTransforms().size()
       *neighbourPatch.patchTransforms().size();

    transforms_.resize(nIntersections);
    patchToPatchList_.resize(nIntersections);

    // Set references to the original patches' point fields
    const pointField& origPatchPoints(ownerPatch.origPatch().localPoints());
    const pointField& nbrOrigPatchPoints
    (
        neighbourPatch.origPatch().localPoints()
    );

    // Create point fields and primitive patches for the transformed patches
    // (originally as a copy of the original patches).
    pointField transfPatchPoints(ownerPatch.origPatch().localPoints());
    pointField transfNbrPatchPoints(neighbourPatch.origPatch().localPoints());

    primitivePatch ownPPatch
    (
        SubList<face>
        (
            ownerPatch.origPatch().localFaces(),
            ownerPatch.origPatch().size()
        ),
        transfPatchPoints
    );
    primitivePatch nbrPPatch
    (
        SubList<face>
        (
            neighbourPatch.origPatch().localFaces(),
            neighbourPatch.origPatch().size()
        ),
        transfNbrPatchPoints
    );

    // Generate patch-to-patch intersections for all possible combinations
    // of rotational couples between original owner and neighbour patches
    label transfIdx = 0;
    forAll(ownerPatch.patchTransforms(), ownTransformi)
    {
        // Transform the original owner patch
        const transformer& ownT = ownerPatch.patchTransforms()[ownTransformi];
        ownT.transformPosition(transfPatchPoints, origPatchPoints);
        ownPPatch.clearGeom();

        forAll(neighbourPatch.patchTransforms(), nbrTransformi)
        {
            // Transform the original neighbour patch
            const transformer& nbrT =
                neighbourPatch.patchTransforms()[nbrTransformi];
            nbrT.transformPosition
            (
                transfNbrPatchPoints,
                nbrOrigPatchPoints
            );
            nbrPPatch.clearGeom();

            // Construct and update patch-to-patch engine for this intersection
            patchToPatchList_.set
            (
                transfIdx,
                new patchToPatches::intersection(false)
            );

            patchToPatchList_[transfIdx].update
            (
                ownPPatch,
                ownPPatch.pointNormals(),
                nbrPPatch,
                transformer::I
            );

            // Append transform pair to the list of transforms
            transforms_[transfIdx] = Pair<transformer>(ownT, nbrT);

            transfIdx++;
        }

        // Reset neighbour patch to the original position for the next rotation
        transfNbrPatchPoints = nbrOrigPatchPoints;
        nbrPPatch.clearGeom();
    }

    // Set source and target coupled faces addressing
    setAddressing();

    // Calculate source and target weights
    calcWeights(true);
    calcWeights(false);

    // Set source and target intersection maps
    setIntersectionMaps();
}


Foam::discreteMixingIntersection::discreteMixingIntersection
(
    const discreteMixingIntersection& fineIntersections,
    const labelList& ownRestrictAddress,
    const labelList& nbrRestrictAddress
)
:
    patchToPatchList_(),
    transforms_(),
    srcAreas_(),
    tgtAreas_(),
    srcAddress_(),
    tgtAddress_(),
    srcWeights_(),
    tgtWeights_(),
    srcZeroWeightFaces_(),
    tgtZeroWeightFaces_(),
    srcIntersectionMap_(),
    tgtIntersectionMap_()
{
    if
    (
        fineIntersections.srcAddress().size() != ownRestrictAddress.size()
     || fineIntersections.tgtAddress().size() != nbrRestrictAddress.size()
    )
    {
        FatalErrorInFunction
            << "Size mismatch between fine interface and agglomeration "
            << "restriction sizes" << nl
            << "Source interface size: "
            << fineIntersections.srcAddress().size() << nl
            << "Source agglomeration restriction size: "
            << ownRestrictAddress.size() << nl
            << "Target interface size: "
            << fineIntersections.tgtAddress().size() << nl
            << "Target agglomeration restriction size: "
            << nbrRestrictAddress.size()
            << exit(FatalError);
    }

    // Agglomerate source areas, addressing, weights and zero-weight faces.
    agglomerate
    (
        fineIntersections.srcAreas(),
        fineIntersections.srcAddress(),
        fineIntersections.srcWeights(),
        ownRestrictAddress,
        nbrRestrictAddress,
        true
    );

    // Agglomerate target areas, addressing, weights and zero-weight faces.
    agglomerate
    (
        fineIntersections.tgtAreas(),
        fineIntersections.tgtAddress(),
        fineIntersections.tgtWeights(),
        nbrRestrictAddress,
        ownRestrictAddress,
        false
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::discreteMixingIntersection::~discreteMixingIntersection()
{}


// ************************************************************************* //
