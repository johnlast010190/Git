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
    (c) 2015 OpenCFD Ltd.
    (c) 2018-2020 OpenFOAM Foundation
    (c) 2021-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "discreteMixingPlanePolyPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "AMIInterpolation/patches/cyclicPeriodicAMI/cyclicPeriodicAMIPolyPatch/cyclicPeriodicAMIPolyPatch.H"
#include "global/constants/mathematical/mathematicalConstants.H"

// For debugging
#include "surfaceFormats/obj/OBJstream.H"
#include "meshes/primitiveMesh/PatchTools/PatchTools.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(discreteMixingPlanePolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, discreteMixingPlanePolyPatch, word);
    addToRunTimeSelectionTable
    (
        polyPatch,
        discreteMixingPlanePolyPatch,
        dictionary
    );


    template<>
    const char* Foam::NamedEnum
    <
        discreteMixingPlanePolyPatch::sectorDefinitions, 2
    >::names[] =
    {
        "periodicPatch",
        "userDefined"
    };

    const NamedEnum
    <
        discreteMixingPlanePolyPatch::sectorDefinitions, 2
    > discreteMixingPlanePolyPatch::sectorDefinitionsTypes_;
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::discreteMixingPlanePolyPatch::resetAMI
(
    const AMIInterpolation::interpolationMethod& AMIMethod
) const
{

    computeSector();

    if (owner())
    {
        srcAMITransforms_.resize(0);
        tgtAMITransforms_.resize(0);

        const transformer oneTr = transformer::I;

        if (sectorDefinition_ == periodicPatch)
        {
            if (!is360())
            {
                // Synchronise the transforms
                syncTransforms();
            }
        }

        // Create copies of both patches' points, transformed to the owner
        pointField thisPoints0(localPoints());
        pointField nbrPoints0(nbrPatch().localPoints());
        transform().transformPosition(nbrPoints0, nbrPoints0);

        autoPtr<OBJstream> ownStr;
        autoPtr<OBJstream> nbrStr;
        if (false)
        {
            const Time& runTime = boundaryMesh().mesh().time();

            fileName dir(runTime.rootPath()/runTime.globalCaseName());
            fileName postfix("_" + runTime.timeName()+"_expanded.obj");

            ownStr.reset(new OBJstream(dir/name() + postfix));
            nbrStr.reset(new OBJstream(dir/nbrPatch().name() + postfix));

            InfoInFunction
                << "patch:" << name()
                << " writing accumulated AMI to " << ownStr().name()
                << " and " << nbrStr().name() << endl;
        }

        // Create another copy
        pointField thisPoints(thisPoints0);
        pointField nbrPoints(nbrPoints0);

        // Create patches for all the points

        // Source patch at initial location
        const primitivePatch thisPatch0
        (
            SubList<face>(localFaces(), size()),
            thisPoints0
        );
        // Source patch that gets moved
        primitivePatch thisPatch
        (
            SubList<face>(localFaces(), size()),
            thisPoints
        );
        // Target patch at initial location
        const primitivePatch nbrPatch0
        (
            SubList<face>(nbrPatch().localFaces(), nbrPatch().size()),
            nbrPoints0
        );
        // Target patch that gets moved
        primitivePatch nbrPatch
        (
            SubList<face>
            (
                this->nbrPatch().localFaces(),
                this->nbrPatch().size()
            ),
            nbrPoints
        );

        {
            label srcSize = returnReduce(thisPatch0.size(), sumOp<label>());
            label tgtSize = returnReduce(nbrPatch0.size(), sumOp<label>());

            Info<< "AMI: Creating addressing and weights between "
                << srcSize << " source faces and " << tgtSize << " target faces"
                << endl;
        }

        // Number of geometry replications
        label iter(0);

        if (ownStr.valid())
        {
            writeOBJ(thisPatch0, ownStr());
        }
        if (nbrStr.valid())
        {
            writeOBJ(nbrPatch0, nbrStr());
        }

        if (debug)
        {
            Info<< "thisSectors(): " << thisSectors() << endl;
            Info<< "neiSectors(): " << nbrSectors() << endl;
        }

        srcAMITransforms_.append(oneTr);
        tgtAMITransforms_.append(oneTr);

        AMIPtr_.reset
        (
            new AMIInterpolation
            (
                thisPatch0,
                nbrPatch0,
                surfPtr(),
                faceAreaIntersect::tmMesh,
                false,
                AMIInterpolation::imPartialFaceAreaWeight,
                AMILowWeightCorrection_,
                AMIReverse_,
                0.01,
                0,
                degToRad
                (
                    debug::floatOptimisationSwitch
                    (
                        "allowedOverlappingAngle",
                        45
                    )
                ),
                1.25,
                false
            )
        );

        // Weight sum averages
        scalar srcSum(0.0);
        scalar tgtSum(0.0);

        if (debug)
        {
            // Weight sum averages
            srcSum =gAverage(AMIPtr_->srcWeightsSum());
            tgtSum = gAverage(AMIPtr_->tgtWeightsSum());
            InfoInFunction
                << "patch:" << name()
                << " srcSum:" << srcSum
                << " tgtSum:" << tgtSum
                << endl;
        }

        nbrPoints = nbrPoints0;
        for (label iNb = 2; iNb <= nbrSectors(); iNb++)
        {
            {
                transformer t = nbrPatchTranformer();
                t.transformPosition(nbrPoints, nbrPoints);
                nbrPatch.clearGeom();
                AMIPtr_->append(thisPatch, nbrPatch, iter+1);

                transformer& tBef = srcAMITransforms_[iter];
                tgtAMITransforms_.append(oneTr);
                srcAMITransforms_.append(t&tBef);
            }

            ++iter;

            if (ownStr.valid())
            {
                writeOBJ(thisPatch, ownStr());
                Info<< "writeOBJon: " << thisPatch.size() << endl;
            }
            if (nbrStr.valid())
            {
                writeOBJ(nbrPatch, nbrStr());
            }


            if (debug)
            {
                srcSum = gAverage(AMIPtr_->srcWeightsSum());
                tgtSum = gAverage(AMIPtr_->tgtWeightsSum());
                InfoInFunction
                    << "patch:" << name()
                    << " iteration:" << iter
                    << " srcSum:" << srcSum
                    << " tgtSum:" << tgtSum
                    << endl;
            }
        }

        for (label iOn = 2; iOn <= thisSectors(); iOn++)
        {
            transformer sF;
            sF = thisPatchTranformer();
            sF.transformPosition(thisPoints, thisPoints);
            thisPatch.clearGeom();
            sF = sF & tgtAMITransforms_[iter];

            for (label iNb = 1; iNb <= nbrSectors(); iNb++)
            {
                //- reset to original position
                if (iNb == 1)
                {
                    nbrPoints = nbrPoints0;
                    nbrPatch.clearGeom();

                    AMIPtr_->append(thisPatch, nbrPatch, iter+1);
                    srcAMITransforms_.append(oneTr);

                }
                else
                {
                    transformer t = nbrPatchTranformer();
                    t.transformPosition(nbrPoints, nbrPoints);
                    nbrPatch.clearGeom();
                    AMIPtr_->append(thisPatch, nbrPatch, iter+1);

                    transformer tF(oneTr);
                    transformer& tBef = srcAMITransforms_[iter];
                    tF = t&tBef;
                    srcAMITransforms_.append(tF);
                }
                tgtAMITransforms_.append(sF);

                ++iter;
                {
                    if (ownStr.valid())
                    {
                        writeOBJ(thisPatch, ownStr());
                    }
                    if (nbrStr.valid())
                    {
                        writeOBJ(nbrPatch, nbrStr());
                    }


                    if (debug)
                    {
                        srcSum = gAverage(AMIPtr_->srcWeightsSum());
                        tgtSum = gAverage(AMIPtr_->tgtWeightsSum());
                        InfoInFunction
                            << "patch:" << name()
                            << " iteration:" << iter
                            << " srcSum:" << srcSum
                            << " tgtSum:" << tgtSum
                            << endl;
                    }
                }
            }
        }


        // Close debug streams
        if (ownStr.valid())
        {
            ownStr.clear();
        }
        if (nbrStr.valid())
        {
            nbrStr.clear();
        }

        // Normalise the weights. Disable printing since weights are
        // still areas.
        AMIPtr_->normaliseWeights(true, false);

        {
            // Print some statistics
            const label nFace = returnReduce
            (
                AMIPtr_->srcWeights().size(), sumOp<label>()
            );

            if (nFace)
            {
                scalarField srcWghtSum(size(), 0);
                forAll(*this, faceI)
                {
                    srcWghtSum[faceI] = sum(AMIPtr_->srcWeights()[faceI]);
                }

                Info<< indent
                    << "AMI: Patch " << name()
                    << " sum(weights) min/max/average = "
                    << gMin(srcWghtSum) << ", "
                    << gMax(srcWghtSum) << ", "
                    << gAverage(srcWghtSum) << endl;
            }
        }
        {
            // Print some statistics
            const label nFace = returnReduce
            (
                AMIPtr_->tgtWeights().size(), sumOp<label>()
            );

            if (nFace)
            {
                scalarField tgtWghtSum(AMIPtr_->tgtWeights().size(), 0);
                forAll(tgtWghtSum, faceI)
                {
                    tgtWghtSum[faceI] = sum(AMIPtr_->tgtWeights()[faceI]);
                }
                Info<< indent
                    << "AMI: Patch " << this->nbrPatch().name()
                    << " sum(weights) min/max/average = "
                    << gMin(tgtWghtSum) << ", "
                    << gMax(tgtWghtSum) << ", "
                    << gAverage(tgtWghtSum) << endl;
            }
        }

        if (false)
        {
            AMIPtr_->visualiseWeights
            (
                *this, this->nbrPatch(),
                false ,
                this->name()+"grgEnd"+this->boundaryMesh().time().timeName()
            );
        }
    }
    if (debug)
    {
        Info<< "srcAMITransforms_ " <<  srcAMITransforms_ <<endl;
        Info<< "tgtAMITransforms_ " <<  tgtAMITransforms_ <<endl;
    }
    AMIPtr_->setTransformations(srcAMITransforms_, tgtAMITransforms_);
}


void Foam::discreteMixingPlanePolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    static_cast<cyclicTransform&>(*this) =
        cyclicTransform
        (
            name(),
            this->primitivePatch::faceAreas(),
            *this,
            nbrPatchName(),
            nbrPatch(),
            matchTolerance()
        );
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::discreteMixingPlanePolyPatch::discreteMixingPlanePolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    cyclicAMIPolyPatch(name, size, start, index, bm, patchType),
    sectorDefinition_(periodicPatch),
    periodicPatchName_(word::null),
    periodicPatchID_(-1),
    sectorMatchTolerance_(matchTolerance()),
    nSectors_(0),
    nSubdivisions_(1),
    srcAMITransforms_(0),
    tgtAMITransforms_(0),
    sectorRotationAxis_(vector::uniform(NaN)),
    sectorRotationCentre_(vector::uniform(NaN))
{
    AMILowWeightCorrection_ = -1.0;
}


Foam::discreteMixingPlanePolyPatch::discreteMixingPlanePolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    cyclicAMIPolyPatch(name, dict, index, bm, patchType),
    sectorDefinition_(periodicPatch),
    periodicPatchName_(dict.lookupOrDefault<word>("periodicPatch", word::null)),
    periodicPatchID_(-1),
    sectorMatchTolerance_
    (
        dict.lookupOrDefault<scalar>("sectorMatchTolerance", matchTolerance())
    ),
    nSectors_(dict.lookupOrDefault<label>("nSectors", 0)),
    nSubdivisions_(dict.lookupOrDefault<label>("nSubdivisions", 1)),
    srcAMITransforms_(0),
    tgtAMITransforms_(0),
    sectorRotationAxis_(vector::uniform(NaN)),
    sectorRotationCentre_(point::uniform(NaN))
{
    AMILowWeightCorrection_ = dict.lookupOrDefault("lowWeightCorrection", 0.01);
    if (dict.found("sectorDefinition"))
    {
        sectorDefinition_ = sectorDefinitionsTypes_.read
        (
            dict.lookup("sectorDefinition")
        );
    }
    if (sectorDefinition_ == userDefined)
    {
        sectorRotationAxis_ = normalised
            (
                dict.lookup<vector>("sectorRotationAxis")
            );
        sectorRotationCentre_ = dict.lookup<point>("sectorRotationCentre");
    }
}


Foam::discreteMixingPlanePolyPatch::discreteMixingPlanePolyPatch
(
    const discreteMixingPlanePolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    cyclicAMIPolyPatch(pp, bm),
    sectorDefinition_(pp.sectorDefinition_),
    periodicPatchName_(pp.periodicPatchName_),
    periodicPatchID_(-1),
    sectorMatchTolerance_(pp.sectorMatchTolerance_),
    nSectors_(pp.nSectors_),
    nSubdivisions_(pp.nSubdivisions_),
    srcAMITransforms_(0),
    tgtAMITransforms_(0),
    sectorRotationAxis_(vector::uniform(NaN)),
    sectorRotationCentre_(vector::uniform(NaN))
{
}


Foam::discreteMixingPlanePolyPatch::discreteMixingPlanePolyPatch
(
    const discreteMixingPlanePolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& nbrPatchName
)
:
    cyclicAMIPolyPatch(pp, bm, index, newSize, newStart, nbrPatchName),
    sectorDefinition_(pp.sectorDefinition_),
    periodicPatchName_(pp.periodicPatchName_),
    periodicPatchID_(-1),
    sectorMatchTolerance_(pp.sectorMatchTolerance_),
    nSectors_(pp.nSectors_),
    nSubdivisions_(pp.nSubdivisions_),
    srcAMITransforms_(0),
    tgtAMITransforms_(0),
    sectorRotationAxis_(vector::uniform(NaN)),
    sectorRotationCentre_(vector::uniform(NaN))
{
}


Foam::discreteMixingPlanePolyPatch::discreteMixingPlanePolyPatch
(
    const discreteMixingPlanePolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    cyclicAMIPolyPatch(pp, bm, index, mapAddressing, newStart),
    sectorDefinition_(pp.sectorDefinition_),
    periodicPatchName_(pp.periodicPatchName_),
    periodicPatchID_(-1),
    sectorMatchTolerance_(pp.sectorMatchTolerance_),
    nSectors_(pp.nSectors_),
    nSubdivisions_(pp.nSubdivisions_),
    srcAMITransforms_(0),
    tgtAMITransforms_(0),
    sectorRotationAxis_(pp.sectorRotationAxis_),
    sectorRotationCentre_(pp.sectorRotationCentre_)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::discreteMixingPlanePolyPatch::write(Ostream& os) const
{
    cyclicAMIPolyPatch::write(os);

    os.writeEntry
    (
        "sectorDefinition", sectorDefinitionsTypes_[sectorDefinition_]
    );

    if (periodicPatchName_ != word::null)
    {
        os.writeEntry("periodicPatch", periodicPatchName_);
    }

    if (sectorMatchTolerance_ != matchTolerance())
    {
        os.writeEntry("sectorMatchTolerance", sectorMatchTolerance_);
    }
    if (nSectors_ != 0)
    {
        os.writeEntry("nSectors", nSectors_);
    }
    os.writeEntry("nSubdivisions", nSubdivisions_);
    if (sectorDefinition_==userDefined)
    {
        os.writeEntry("sectorRotationAxis", sectorRotationAxis_);
        os.writeEntry("sectorRotationCentre", sectorRotationCentre_);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AMIInterpolation/patches/discreteMixingPlane/discreteMixingPlanePolyPatch/discreteMixingPlanePolyPatchSectors.C"

// ************************************************************************* //
