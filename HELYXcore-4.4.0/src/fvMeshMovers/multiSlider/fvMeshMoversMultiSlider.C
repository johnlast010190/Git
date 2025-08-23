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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2017-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "multiSlider/fvMeshMoversMultiSlider.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshMovers
{
    defineTypeNameAndDebug(multiSlider, 0);
    addToRunTimeSelectionTable(fvMeshMover, multiSlider, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshMovers::multiSlider::multiSlider(fvMesh& mesh)
:
    fvMeshMover(mesh),
    coeffs_(dict().subDict(typeName + "Coeffs")),
    displacement_(mesh.points().size(), vector::zero)
{
    zoneIDs_.setSize(coeffs_.size());
    pointIDs_.setSize(coeffs_.size());
    movingIDs_.setSize(coeffs_.size());
    motionProfile_.setSize(coeffs_.size());
    velocityProfile_.setSize(coeffs_.size());
    initialDisplacement_.setSize(coeffs_.size());
    bufferSize_.setSize(coeffs_.size());
    blendType_.setSize(coeffs_.size());

    label zoneI = 0;

    const Time& t = mesh.time();
    scalar time = t.timeOutputValue();

    forAllConstIter(dictionary, coeffs_, iter)
    {
        if (iter().isDict())
        {
            zoneIDs_[zoneI] = mesh.cellZones().findZoneID(iter().keyword());

            if (zoneIDs_[zoneI] == -1)
            {
                FatalIOErrorIn
                (
                    "multiSliderMotionFvMesh::"
                    "multiSliderMotionFvMesh(const IOobject&)",
                    coeffs_
                )   << "Cannot find cellZone named " << iter().keyword()
                    << ". Valid zones are " << mesh.cellZones().names()
                    << exit(FatalIOError);
            }

            const dictionary& subDict = iter().dict();

            // List of moving patch ID's to define fixed mesh
            // displacement region
            wordList patchNames = wordList(subDict.lookup("movingPatches"));
            HashSet<word> patchNameSet(patchNames);

            DynamicList<label> mIDs(patchNames.size());
            forAll(mesh.boundary(), patchI)
            {
                word name = mesh.boundary()[patchI].name();
                if (patchNameSet.found(name))
                {
                    mIDs.append(patchI);
                }
            }
            mIDs.shrink();
            movingIDs_[zoneI].transfer(mIDs);

            // Velocity data entry for each zone
            if (subDict.found("velocity"))
            {
                velocityProfile_[zoneI] = true;
                motionProfile_.set
                (
                    zoneI,
                    Function1<vector>::New("velocity", subDict)
                );
            }
            else if (subDict.found("displacement"))
            {
                velocityProfile_[zoneI] = false;
                motionProfile_.set
                (
                    zoneI,
                    Function1<vector>::New("displacement", subDict)
                );
            }
            else
            {
                FatalErrorInFunction
                    << "Need displacement or velocity profile for zone: "
                    << iter().keyword()
                    << exit(FatalError);
            }

            if (time < SMALL)
            {
                // Define initial mesh displacement if present
                initialDisplacement_[zoneI] =
                    subDict.lookupOrDefault<vector>
                    (
                        "initialDisplacement",
                        vector::zero
                    );
            }
            else
            {
                initialDisplacement_[zoneI] = vector::zero;
            }

            // Define buffer region size and blending function
            bufferSize_[zoneI] =
                subDict.lookupOrDefault<scalar>("bufferSize", 0.01);
            blendType_[zoneI] =
                subDict.lookupOrDefault<word>("blendType", "linear");

            // Collect points of cell zone
            const cellZone& cz = mesh.cellZones()[zoneIDs_[zoneI]];

            boolList movePts(mesh.nPoints(), false);

            forAll(cz, i)
            {
                label cellI = cz[i];
                const cell& c = mesh.cells()[cellI];
                forAll(c, j)
                {
                    const face& f = mesh.faces()[c[j]];
                    forAll(f, k)
                    {
                        label pointI = f[k];
                        movePts[pointI] = true;
                    }
                }
            }

            syncTools::syncPointList(mesh, movePts, orEqOp<bool>(), false);

            DynamicList<label> ptIDs(mesh.nPoints());
            forAll(movePts, i)
            {
                if (movePts[i])
                {
                    ptIDs.append(i);
                }
            }

            pointIDs_[zoneI].transfer(ptIDs);

            Info<< "Applying solid body motion to "
                << pointIDs_[zoneI].size() << " points of cellZone "
                << iter().keyword() << endl;

            zoneI++;
        }
    }

    zoneIDs_.setSize(zoneI);
    pointIDs_.setSize(zoneI);
    movingIDs_.setSize(zoneI);
    motionProfile_.setSize(zoneI);
    velocityProfile_.setSize(zoneI);
    initialDisplacement_.setSize(zoneI);
    bufferSize_.setSize(zoneI);
    blendType_.setSize(zoneI);

    // Move points according if initial displacement defined
    slidePoints(initialDisplacement_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshMovers::multiSlider::~multiSlider()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::fvMeshMovers::multiSlider::makePatch(const labelList& patchIDs)
{
    const polyBoundaryMesh& patches = mesh().boundaryMesh();

    // Count faces
    label nFaces = 0;

    forAll(patchIDs, i)
    {
        const polyPatch& pp = patches[patchIDs[i]];

        nFaces += pp.size();
    }

    // Collect faces
    labelList addressing(nFaces);
    nFaces = 0;

    forAll(patchIDs, i)
    {
        const polyPatch& pp = patches[patchIDs[i]];

        label meshFaceI = pp.start();

        forAll(pp, i)
        {
            addressing[nFaces++] = meshFaceI++;
        }
    }

    return autoPtr<indirectPrimitivePatch>
    (
        new indirectPrimitivePatch
        (
            IndirectList<face>(mesh().faces(), addressing),
            mesh().points()
        )
    );
}


void Foam::fvMeshMovers::multiSlider::slidePoints
(
    const vectorField& zoneDisplacements
)
{
    forAll(zoneIDs_, i)
    {
        const labelList& zonePoints = pointIDs_[i];
        if (mag(zoneDisplacements[i]) > SMALL)
        {
            vector dir = zoneDisplacements[i]/mag(zoneDisplacements[i]);

            const labelList& movingPatchIds = movingIDs_[i];
            autoPtr<indirectPrimitivePatch> pp(makePatch(movingPatchIds));

            scalar minRegionDP = GREAT;
            scalar maxRegionDP = -GREAT;

            forAll(pp().localPoints(), pointI)
            {
                scalar dProd = (pp().localPoints()[pointI] & dir);
                minRegionDP =  min(minRegionDP, dProd);
                maxRegionDP =  max(maxRegionDP, dProd);
            }

            scalar minZoneDP = GREAT;
            scalar maxZoneDP = -GREAT;

            pointField zPts(mesh().points(), zonePoints);
            scalarField zDP(zPts.size(), 0.0);

            forAll(zPts, pointI)
            {
                zDP[pointI] = (zPts[pointI] & dir);
                minZoneDP = min(minZoneDP, zDP[pointI]);
                maxZoneDP = max(maxZoneDP, zDP[pointI]);
            }

            Foam::reduce
            (
                std::tie(minRegionDP, maxRegionDP, minZoneDP, maxZoneDP),
                ParallelOp
                <
                    minOp<scalar>, maxOp<scalar>, minOp<scalar>, maxOp<scalar>
                >{},
                mesh().comm()
            );

            scalarField dispRatio(zPts.size(), 0.0);

            scalar totalLength = maxZoneDP - minZoneDP;
            scalar bufferLength = bufferSize_[i] * totalLength;

            if
            (
                bufferLength > 0.25*(maxZoneDP - maxRegionDP)
             || bufferLength > 0.25*(minRegionDP - minZoneDP)
            )
            {
                bufferLength = 0.0;
            }

            forAll(zPts, pointI)
            {
                scalar dn = zDP[pointI];

                if
                (
                    dn > minZoneDP + bufferLength
                 && dn < maxZoneDP - bufferLength
                )
                {
                    if
                    (
                        dn >= minRegionDP - bufferLength
                     && dn <= maxRegionDP + bufferLength
                    )
                    {
                        dispRatio[pointI] = 1.0;
                    }
                    else
                    {
                        if (dn > maxRegionDP + bufferLength)
                        {
                            scalar num = dn - (maxRegionDP + bufferLength);
                            scalar den =
                                maxZoneDP - maxRegionDP - 2*bufferLength;
                            dispRatio[pointI] = 1.0 - num/den;
                        }
                        else
                        {
                            scalar num = dn - (minZoneDP + bufferLength);
                            scalar den =
                                minRegionDP - minZoneDP - 2*bufferLength;
                            dispRatio[pointI] = num/den;
                        }
                    }
                }
            }

            if (blendType_[i] == "linear")
            {
                UIndirectList<point>(displacement_, zonePoints) =
                    dispRatio*zoneDisplacements[i];
            }
            else if (blendType_[i] == "cos")
            {
                UIndirectList<point>(displacement_, zonePoints) =
                    0.5*
                    (
                        1.0 - Foam::cos
                        (
                            dispRatio*Foam::constant::mathematical::pi
                        )
                    )*zoneDisplacements[i];
            }
            else if (blendType_[i] == "smoothStep")
            {
                UIndirectList<point>(displacement_, zonePoints) =
                    (3.0*pow(dispRatio, 2) - 2.0*pow(dispRatio, 3))
                   *zoneDisplacements[i];
            }
            else
            {
                FatalErrorInFunction
                    << "Blending scheme " << blendType_[i] << " not available."
                    << exit(FatalError);
            }
        }
        else
        {
            UIndirectList<point>(displacement_, zonePoints) = vector::zero;
        }
    }

    syncTools::syncPointList
    (
        mesh(),
        displacement_,
        minMagSqrEqOp<point>(),         // combine op
        vector(GREAT, GREAT, GREAT)     // null value
    );

    mesh().movePoints(mesh().points() + displacement_);
}


bool Foam::fvMeshMovers::multiSlider::update()
{
    static bool hasWarned = false;

    const Time& t = mesh().time();
    vectorField zoneDisplacements(zoneIDs_.size());
    forAll(zoneIDs_, i)
    {
        scalar t0 = t.timeOutputValue() - t.deltaTValue();
        if (velocityProfile_[i])
        {
            vector vel = motionProfile_[i].value(t0);
            zoneDisplacements[i] = vel*t.deltaTValue();
        }
        else
        {
            scalar t1 = t.timeOutputValue();
            vector disp0 = motionProfile_[i].value(t0);
            vector disp1 = motionProfile_[i].value(t1);
            zoneDisplacements[i] = disp1 - disp0;
        }
    }

    slidePoints(zoneDisplacements);

    if (mesh().foundObject<volVectorField>("U"))
    {
        const_cast<volVectorField&>(mesh().lookupObject<volVectorField>("U"))
            .correctBoundaryConditions();
    }
    else if (!hasWarned)
    {
        hasWarned = true;

        WarningInFunction
            << "Did not find volVectorField U."
            << " Not updating U boundary conditions." << endl;
    }

    return true;
}


// ************************************************************************* //
