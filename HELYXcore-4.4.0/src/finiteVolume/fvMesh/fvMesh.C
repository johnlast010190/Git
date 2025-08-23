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
    (c) 2016-2017 OpenCFD Ltd.
    (c) 2011-2024 OpenFOAM Foundation
    (c) 2010-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/GeometricFields/pointFields/pointFields.H"
#include "fields/cloud/mapClouds.H"
#include "fields/fvsPatchFields/derived/polyFaces/polyFacesFvsPatchLabelField.H"
#include "fields/fvsPatchFields/derived/nonConformalPolyFaces/nonConformalPolyFacesFvsPatchLabelField.H"
#include "interpolation/mapping/fvFieldMappers/MapFvFields.H"
#include "meshes/pointMesh/pointMeshMapper/pointMeshMapper.H"
#include "meshes/pointMesh/pointMeshMapper/MapPointField.H"
#include "meshes/MeshObject/MeshObject.H"
#include "fvMesh/fvMeshLduAddressing.H"
#include "fvMesh/fvMeshMapper/fvMeshMapper.H"
#include "fvMesh/fvMeshStitchers/fvMeshStitcher/fvMeshStitcher.H"
#include "fvMesh/fvMeshMovers/fvMeshMover/fvMeshMover.H"
#include "fvMesh/fvMeshTopoChangers/fvMeshTopoChanger/fvMeshTopoChanger.H"
#include "fvMesh/fvMeshDistributors/fvMeshDistributor/fvMeshDistributor.H"
#include "fvMesh/fvMeshGIBChangers/fvMeshGIBChanger/fvMeshGIBChanger.H"
#include "fvMesh/fvPatches/derived/nonConformalOrig/nonConformalOrigFvPatch.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "VectorN/finiteVolume/fields/volFields/volVectorNFields.H"
#include "nonConformal/nonConformalFuncs/nonConformalFuncs.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMesh, 0);
}

const Foam::HashSet<Foam::word> Foam::fvMesh::geometryFields
{
    "Vc",
    "Vc0",
    "Vc00",
    "Sf",
    "magSf",
    "Cc",
    "Cf",
    "meshPhi"
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMesh::clearGeomNotOldVol()
{
    if (debug)
    {
        InfoInFunction << "clearGeomNotOldVol" << endl;
    }

    meshObject::clearUpto
    <
        fvMesh,
        GeometricMeshObject,
        MoveableMeshObject
    >(*this);

    meshObject::clearUpto
    <
        lduMesh,
        GeometricMeshObject,
        MoveableMeshObject
    >(*this);

    deleteDemandDrivenData(VPtr_);
    deleteDemandDrivenData(SfSlicePtr_);
    deleteDemandDrivenData(SfPtr_);
    deleteDemandDrivenData(magSfSlicePtr_);
    deleteDemandDrivenData(magSfPtr_);
    deleteDemandDrivenData(CSlicePtr_);
    deleteDemandDrivenData(CPtr_);
    deleteDemandDrivenData(CfSlicePtr_);
    deleteDemandDrivenData(CfPtr_);
}


void Foam::fvMesh::updateGeomNotOldVol()
{
    bool haveV = (VPtr_ != nullptr);
    bool haveSf = (SfSlicePtr_ != nullptr || SfPtr_ != nullptr);
    bool haveMagSf = (magSfSlicePtr_ != nullptr || magSfPtr_ != nullptr);
    bool haveCP = (CSlicePtr_ != nullptr || CPtr_ != nullptr);
    bool haveCf = (CfSlicePtr_ != nullptr || CfPtr_ != nullptr);

    clearGeomNotOldVol();

    // Now recreate the fields
    if (haveV)
    {
        (void)V();
    }
    if (haveSf)
    {
        (void)Sf();
    }
    if (haveMagSf)
    {
        (void)magSf();
    }
    if (haveCP)
    {
        (void)C();
    }
    if (haveCf)
    {
        (void)Cf();
    }
}


void Foam::fvMesh::clearGeom()
{
    if (debug)
    {
        InfoInFunction << "Clearing geometric data" << endl;
    }

    clearGeomNotOldVol();

    deleteDemandDrivenData(phiPtr_);
    deleteDemandDrivenData(V0Ptr_);
    deleteDemandDrivenData(V00Ptr_);
}


void Foam::fvMesh::clearAddressing(const bool isMeshUpdate)
{
    if (debug)
    {
        InfoInFunction << "isMeshUpdate: " << isMeshUpdate << endl;
    }

    if (isMeshUpdate)
    {
        // Part of a mesh update. Keep meshObjects that have an topoChange
        // callback
        meshObject::clearUpto
        <
            fvMesh,
            TopologicalMeshObject,
            UpdateableMeshObject
        >
        (
            *this
        );
        meshObject::clearUpto
        <
            lduMesh,
            TopologicalMeshObject,
            UpdateableMeshObject
        >
        (
            *this
        );
    }
    else
    {
        meshObject::clear<fvMesh, TopologicalMeshObject>(*this);
        meshObject::clear<lduMesh, TopologicalMeshObject>(*this);
    }

    deleteDemandDrivenData(lduPtr_);
    deleteDemandDrivenData(polyFacesBfPtr_);
    deleteDemandDrivenData(polyBFaceOffsetsPtr_);
    deleteDemandDrivenData(polyBFaceOffsetPatchesPtr_);
    deleteDemandDrivenData(polyBFaceOffsetPatchFacesPtr_);
    deleteDemandDrivenData(polyBFacePatchesPtr_);
    deleteDemandDrivenData(polyBFacePatchFacesPtr_);
    deleteDemandDrivenData(ownerBfPtr_);
}


void Foam::fvMesh::storeOldTimeFields()
{
    storeOldTimeFields<PointField>();
    storeOldTimeFields<VolField>();
    storeOldTimeFields<SurfaceField>();
}


void Foam::fvMesh::nullOldestTimeFields()
{
    nullOldestTimeFields<PointField>();
    nullOldestTimeFields<VolField>();
    nullOldestTimeFields<SurfaceField>();
}


Foam::wordList Foam::fvMesh::polyFacesPatchTypes() const
{
    wordList wantedPatchTypes
    (
        boundary().size(),
        polyFacesFvsPatchLabelField::typeName
    );

    forAll(boundary(), patchi)
    {
        const fvPatch& fvp = boundary()[patchi];

        if (isA<nonConformalFvPatch>(fvp))
        {
            wantedPatchTypes[patchi] =
                nonConformalPolyFacesFvsPatchLabelField::typeName;
        }
    }

    return wantedPatchTypes;
}


Foam::surfaceLabelField::Boundary& Foam::fvMesh::polyFacesBfRef()
{
    if (!polyFacesBfPtr_)
    {
        polyFacesBf();
    }

    return *polyFacesBfPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMesh::fvMesh
(
    const IOobject& io,
    const bool changers,
    const stitchType stitch
)
:
    polyMesh(io),
    surfaceInterpolation(*this),
    data(static_cast<const objectRegistry&>(*this)),
    boundary_(*this, boundaryMesh()),
    stitcher_(fvMeshStitcher::New(*this, changers)),
    topoChanger_(nullptr),
    distributor_(nullptr),
    mover_(nullptr),
    GIBChanger_(nullptr),
    lduPtr_(nullptr),
    polyFacesBfPtr_(nullptr),
    polyBFaceOffsetsPtr_(nullptr),
    polyBFaceOffsetPatchesPtr_(nullptr),
    polyBFaceOffsetPatchFacesPtr_(nullptr),
    polyBFacePatchesPtr_(nullptr),
    polyBFacePatchFacesPtr_(nullptr),
    ownerBfPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    VPtr_(nullptr),
    V0Ptr_(nullptr),
    V00Ptr_(nullptr),
    SfSlicePtr_(nullptr),
    SfPtr_(nullptr),
    magSfSlicePtr_(nullptr),
    magSfPtr_(nullptr),
    CSlicePtr_(nullptr),
    CPtr_(nullptr),
    CfSlicePtr_(nullptr),
    CfPtr_(nullptr),
    phiPtr_(nullptr)
{
    if (debug)
    {
        InfoInFunction << "Constructing fvMesh from IOobject" << endl;
    }

    if (changers && needsStitching())
    {
        // Construct non-conformal couples on-the-fly
        nonConformalFuncs::createNonConformalCouples(*this);
    }
    else if (stitch != stitchType::none)
    {
        // Stitch or re-stitch if necessary
        stitcher_->connect(false, stitch == stitchType::geometric, true);
    }

    // Construct changers
    if (changers)
    {
        topoChanger_.set(fvMeshTopoChanger::New(*this).ptr());
        distributor_.set(fvMeshDistributor::New(*this).ptr());
        mover_.set(fvMeshMover::New(*this).ptr());
        GIBChanger_.set(fvMeshGIBChanger::New(*this).ptr());
    }

    // Check the existence of the cell volumes and read if present
    if (fileHandler().isFile(time().timePath()/"Vc0"))
    {
        V0Ptr_ = new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "Vc0",
                time().timeName(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                true
            ),
            *this
        );
    }

    // Check the existence of the mesh fluxes and read if present
    if (fileHandler().isFile(time().timePath()/"meshPhi"))
    {
        phiPtr_ = new surfaceScalarField
        (
            IOobject
            (
                "meshPhi",
                time().timeName(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                true
            ),
            *this,
            true,
            false
        );
    }
}


Foam::fvMesh::fvMesh
(
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    labelList&& allOwner,
    labelList&& allNeighbour,
    const bool syncPar
)
:
    polyMesh
    (
        io,
        std::move(points),
        std::move(faces),
        std::move(allOwner),
        std::move(allNeighbour),
        syncPar
    ),
    surfaceInterpolation(*this),
    data(static_cast<const objectRegistry&>(*this)),
    boundary_(*this, boundaryMesh()),
    stitcher_(nullptr),
    topoChanger_(nullptr),
    distributor_(nullptr),
    mover_(nullptr),
    GIBChanger_(nullptr),
    lduPtr_(nullptr),
    polyFacesBfPtr_(nullptr),
    polyBFaceOffsetsPtr_(nullptr),
    polyBFaceOffsetPatchesPtr_(nullptr),
    polyBFaceOffsetPatchFacesPtr_(nullptr),
    polyBFacePatchesPtr_(nullptr),
    polyBFacePatchFacesPtr_(nullptr),
    ownerBfPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    VPtr_(nullptr),
    V0Ptr_(nullptr),
    V00Ptr_(nullptr),
    SfSlicePtr_(nullptr),
    SfPtr_(nullptr),
    magSfSlicePtr_(nullptr),
    magSfPtr_(nullptr),
    CSlicePtr_(nullptr),
    CPtr_(nullptr),
    CfSlicePtr_(nullptr),
    CfPtr_(nullptr),
    phiPtr_(nullptr)
{
    if (debug)
    {
        InfoInFunction << "Constructing fvMesh from components" << endl;
    }
}


Foam::fvMesh::fvMesh
(
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    cellList&& cells,
    const bool syncPar,
    const bool autoWrite
)
:
    polyMesh
    (
        io,
        std::move(points),
        std::move(faces),
        std::move(cells),
        syncPar,
        autoWrite
    ),
    surfaceInterpolation(*this),
    data(static_cast<const objectRegistry&>(*this)),
    boundary_(*this),
    stitcher_(nullptr),
    topoChanger_(nullptr),
    distributor_(nullptr),
    mover_(nullptr),
    GIBChanger_(nullptr),
    lduPtr_(nullptr),
    polyFacesBfPtr_(nullptr),
    polyBFaceOffsetsPtr_(nullptr),
    polyBFaceOffsetPatchesPtr_(nullptr),
    polyBFaceOffsetPatchFacesPtr_(nullptr),
    polyBFacePatchesPtr_(nullptr),
    polyBFacePatchFacesPtr_(nullptr),
    ownerBfPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    VPtr_(nullptr),
    V0Ptr_(nullptr),
    V00Ptr_(nullptr),
    SfSlicePtr_(nullptr),
    SfPtr_(nullptr),
    magSfSlicePtr_(nullptr),
    magSfPtr_(nullptr),
    CSlicePtr_(nullptr),
    CPtr_(nullptr),
    CfSlicePtr_(nullptr),
    CfPtr_(nullptr),
    phiPtr_(nullptr)
{
    if (debug)
    {
        InfoInFunction << "Constructing fvMesh from components" << endl;
    }
}


Foam::fvMesh::fvMesh
(
    const IOobject& io,
    pointField&& points,
    const cellShapeList& shapes,
    const faceListList& boundaryFaces,
    const wordList& boundaryPatchNames,
    const PtrList<dictionary>& boundaryDicts,
    const word& defaultBoundaryPatchName,
    const word& defaultBoundaryPatchType,
    const bool syncPar
)
:
    polyMesh
    (
       io,
       std::move(points),
       shapes,
       boundaryFaces,
       boundaryPatchNames,
       boundaryDicts,
       defaultBoundaryPatchName,
       defaultBoundaryPatchType,
       syncPar
    ),
    surfaceInterpolation(*this),
    data(static_cast<const objectRegistry&>(*this)),
    boundary_(*this,boundaryMesh()),
    stitcher_(nullptr),
    topoChanger_(nullptr),
    distributor_(nullptr),
    mover_(nullptr),
    GIBChanger_(nullptr),
    lduPtr_(nullptr),
    polyFacesBfPtr_(nullptr),
    polyBFaceOffsetsPtr_(nullptr),
    polyBFaceOffsetPatchesPtr_(nullptr),
    polyBFaceOffsetPatchFacesPtr_(nullptr),
    polyBFacePatchesPtr_(nullptr),
    polyBFacePatchFacesPtr_(nullptr),
    ownerBfPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    VPtr_(nullptr),
    V0Ptr_(nullptr),
    V00Ptr_(nullptr),
    SfSlicePtr_(nullptr),
    SfPtr_(nullptr),
    magSfSlicePtr_(nullptr),
    magSfPtr_(nullptr),
    CSlicePtr_(nullptr),
    CPtr_(nullptr),
    CfSlicePtr_(nullptr),
    CfPtr_(nullptr),
    phiPtr_(nullptr)
{
    if (debug)
    {
        InfoInFunction << "Constructing fvMesh from blockMesh" << endl;
    }
}


Foam::fvMesh::fvMesh(fvMesh&& mesh)
:
    polyMesh(std::move(mesh)),
    surfaceInterpolation(std::move(mesh)),
    data(static_cast<const objectRegistry&>(*this)),
    boundary_(std::move(mesh.boundary_)),
    stitcher_(std::move(mesh.stitcher_)),
    topoChanger_(std::move(mesh.topoChanger_)),
    distributor_(std::move(mesh.distributor_)),
    mover_(std::move(mesh.mover_)),
    GIBChanger_(std::move(mesh.GIBChanger_)),
    lduPtr_(std::move(mesh.lduPtr_)),
    polyFacesBfPtr_(std::move(mesh.polyFacesBfPtr_)),
    polyBFaceOffsetsPtr_(std::move(mesh.polyBFaceOffsetsPtr_)),
    polyBFaceOffsetPatchesPtr_(std::move(mesh.polyBFaceOffsetPatchesPtr_)),
    polyBFaceOffsetPatchFacesPtr_
    (
        std::move(mesh.polyBFaceOffsetPatchFacesPtr_)
    ),
    polyBFacePatchesPtr_(std::move(mesh.polyBFacePatchesPtr_)),
    polyBFacePatchFacesPtr_(std::move(mesh.polyBFacePatchFacesPtr_)),
    ownerBfPtr_(std::move(mesh.ownerBfPtr_)),
    curTimeIndex_(mesh.curTimeIndex_),
    VPtr_(std::move(mesh.VPtr_)),
    V0Ptr_(std::move(mesh.V0Ptr_)),
    V00Ptr_(std::move(mesh.V00Ptr_)),
    SfSlicePtr_(std::move(mesh.SfSlicePtr_)),
    SfPtr_(std::move(mesh.SfPtr_)),
    magSfSlicePtr_(std::move(mesh.magSfSlicePtr_)),
    magSfPtr_(std::move(mesh.magSfPtr_)),
    CSlicePtr_(std::move(mesh.CSlicePtr_)),
    CPtr_(std::move(mesh.CPtr_)),
    CfSlicePtr_(std::move(mesh.CfSlicePtr_)),
    CfPtr_(std::move(mesh.CfPtr_)),
    phiPtr_(std::move(mesh.phiPtr_))
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << "Moving fvMesh" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMesh::~fvMesh()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMesh::addFvPatches
(
    const List<polyPatch*> & p,
    const bool validBoundary
)
{
    if (boundary().size())
    {
        FatalErrorInFunction
            << " boundary already exists"
            << abort(FatalError);
    }

    // first add polyPatches
    addPatches(p, validBoundary);
    boundary_.addPatches(boundaryMesh());
}


Foam::polyMesh::readUpdateState Foam::fvMesh::readUpdate
(
    const stitchType stitch
)
{
    if (debug)
    {
        InfoInFunction << "Updating fvMesh.  ";
    }

    polyMesh::readUpdateState state = polyMesh::readUpdate();

    if (state == polyMesh::TOPO_PATCH_CHANGE)
    {
        boundary_.readUpdate(boundaryMesh());
    }

    if
    (
        stitcher_.valid()
     && stitch != stitchType::none
     && state != polyMesh::UNCHANGED
    )
    {
        stitcher_->disconnect(false, stitch == stitchType::geometric);
    }

    if (state == polyMesh::TOPO_PATCH_CHANGE)
    {
        if (debug)
        {
            Info<< "Boundary and topological update" << endl;
        }

        clearOut();

    }
    else if (state == polyMesh::TOPO_CHANGE)
    {
        if (debug)
        {
            Info<< "Topological update" << endl;
        }

        clearOut();
    }
    else if (state == polyMesh::POINTS_MOVED)
    {
        if (debug)
        {
            Info<< "Point motion update" << endl;
        }

        clearGeom();
    }
    else
    {
        if (debug)
        {
            Info<< "No update" << endl;
        }
    }

    if
    (
        stitcher_.valid()
     && stitch != stitchType::none
     && state != polyMesh::UNCHANGED
    )
    {
        stitcher_->connect(false, stitch == stitchType::geometric, true);
    }

    // If the mesh has been re-stitched with different geometry, then the
    // finite-volume topology has changed
    if
    (
        stitcher_.valid()
     && stitcher_->stitches()
     && state == polyMesh::POINTS_MOVED
    )
    {
        state = polyMesh::TOPO_CHANGE;
    }

    return state;
}


const Foam::lduAddressing& Foam::fvMesh::lduAddr() const
{
    if (!lduPtr_)
    {
        lduPtr_ = new fvMeshLduAddressing(*this);
    }

    return *lduPtr_;
}


Foam::SolverPerformance<Foam::scalar> Foam::fvMesh::solve
(
    fvMatrix<scalar>& m,
    const dictionary& dict
) const
{
    // Redirect to fvMatrix solver
    return m.solveSegregatedOrCoupled(dict);
}


Foam::SolverPerformance<Foam::vector> Foam::fvMesh::solve
(
    fvMatrix<vector>& m,
    const dictionary& dict
) const
{
    // Redirect to fvMatrix solver
    return m.solveSegregatedOrCoupled(dict);
}


Foam::SolverPerformance<Foam::sphericalTensor> Foam::fvMesh::solve
(
    fvMatrix<sphericalTensor>& m,
    const dictionary& dict
) const
{
    // Redirect to fvMatrix solver
    return m.solveSegregatedOrCoupled(dict);
}


Foam::SolverPerformance<Foam::symmTensor> Foam::fvMesh::solve
(
    fvMatrix<symmTensor>& m,
    const dictionary& dict
) const
{
    // Redirect to fvMatrix solver
    return m.solveSegregatedOrCoupled(dict);
}


Foam::SolverPerformance<Foam::tensor> Foam::fvMesh::solve
(
    fvMatrix<tensor>& m,
    const dictionary& dict
) const
{
    // Redirect to fvMatrix solver
    return m.solveSegregatedOrCoupled(dict);
}


Foam::IOobject Foam::fvMesh::polyFacesBfIO(const IOobject::readOption r) const
{
    return
        IOobject
        (
            "polyFaces",
            pointsInstance(),
            typeName,
            *this,
            r,
            IOobject::NO_WRITE,
            false
        );
}


const Foam::surfaceLabelField::Boundary& Foam::fvMesh::polyFacesBf() const
{
    if (!polyFacesBfPtr_)
    {
        polyFacesBfPtr_ =
            new surfaceLabelField::Boundary
            (
                boundary(),
                surfaceLabelField::null(),
                polyFacesPatchTypes(),
                boundaryMesh().types()
            );
    }

    return *polyFacesBfPtr_;
}


const Foam::CompactListList<Foam::label>&
Foam::fvMesh::polyBFacePatches() const
{
    if (!polyBFacePatchesPtr_)
    {
        const label nPolyBFaces = nFaces() - nInternalFaces();

        // Count face-poly-bFaces to get the offsets
        polyBFaceOffsetsPtr_ = new labelList(nPolyBFaces + 1, 0);
        labelList& offsets = *polyBFaceOffsetsPtr_;
        forAll(boundary(), patchi)
        {
            if (!isA<indirectPolyPatch>(boundary()[patchi].patch()))
            {
                forAll(boundary()[patchi], patchFacei)
                {
                    const label polyBFacei =
                        (
                            polyFacesBfPtr_
                        ? (*polyFacesBfPtr_)[patchi][patchFacei]
                        : boundary()[patchi].start() + patchFacei
                        )
                    - nInternalFaces();
                    offsets[polyBFacei + 1] ++;
                }
            }
        }
        for (label polyBFacei = 0; polyBFacei < nPolyBFaces; ++ polyBFacei)
        {
            offsets[polyBFacei + 1] += offsets[polyBFacei];
        }

        // Set the poly-bFace patches and patch-faces, using the offsets as
        // counters
        polyBFaceOffsetPatchesPtr_ = new labelList(offsets.last());
        polyBFaceOffsetPatchFacesPtr_ = new labelList(offsets.last());
        labelUList& patches = *polyBFaceOffsetPatchesPtr_;
        labelUList& patchFaces = *polyBFaceOffsetPatchFacesPtr_;
        forAll(boundary(), patchi)
        {
            if (!isA<indirectPolyPatch>(boundary()[patchi].patch()))
            {
                forAll(boundary()[patchi], patchFacei)
                {
                    const label polyBFacei =
                        (
                            polyFacesBfPtr_
                            ? (*polyFacesBfPtr_)[patchi][patchFacei]
                            : boundary()[patchi].start() + patchFacei
                        )
                        - nInternalFaces();
                    patches[offsets[polyBFacei]] = patchi;
                    patchFaces[offsets[polyBFacei]] = patchFacei;
                    offsets[polyBFacei] ++;
                }
            }
        }

        // Restore the offsets by removing the count
        for
        (
            label polyBFacei = nPolyBFaces - 1;
            polyBFacei >= 0;
            -- polyBFacei
        )
        {
            offsets[polyBFacei + 1] = offsets[polyBFacei];
        }
        offsets[0] = 0;

        // List-lists
        polyBFacePatchesPtr_ =
            new CompactListList<label>(offsets, patches);
        polyBFacePatchFacesPtr_ =
            new CompactListList<label>(offsets, patchFaces);
    }

    return *polyBFacePatchesPtr_;
}


const Foam::CompactListList<Foam::label>&
Foam::fvMesh::polyBFacePatchFaces() const
{
    if (!polyBFacePatchFacesPtr_)
    {
        polyBFacePatches();
    }

    return *polyBFacePatchFacesPtr_;
}


const Foam::surfaceLabelField::Boundary& Foam::fvMesh::ownerBf() const
{
    if (!ownerBfPtr_)
    {
        ownerBfPtr_ =
            new surfaceLabelField::Boundary
            (
                boundary(),
                surfaceLabelField::null(),
                wordList
                (
                    boundary().size(),
                    calculatedFvsPatchLabelField::typeName
                ),
                boundaryMesh().types()
            );

        forAll(boundary(), patchi)
        {
            (*ownerBfPtr_)[patchi] =
                labelField(faceOwner(), polyFacesBf()[patchi]);
        }
    }

    return *ownerBfPtr_;
}


const Foam::fvMeshStitcher& Foam::fvMesh::stitcher() const
{
    return stitcher_();
}


const Foam::fvMeshTopoChanger& Foam::fvMesh::topoChanger() const
{
    return topoChanger_();
}


const Foam::fvMeshDistributor& Foam::fvMesh::distributor() const
{
    return distributor_();
}


const Foam::fvMeshMover& Foam::fvMesh::mover() const
{
    return mover_();
}


const Foam::fvMeshGIBChanger& Foam::fvMesh::GIBChanger() const
{
    return GIBChanger_();
}


bool Foam::fvMesh::hasChangers() const
{
    return
        topoChanger_.valid()
     || mover_.valid()
     || distributor_.valid()
     || GIBChanger_.valid();
}


bool Foam::fvMesh::needsStitching() const
{
    bool hasNonConformalCouples = false;

    forAll(boundary(), patchi)
    {
        if (isA<nonConformalOrigFvPatch>(boundary()[patchi]))
        {
            hasNonConformalCouples = true;
            break;
        }
    }

    IOobject nccDictIO
    (
        "nonConformalCouplesDict",
        this->time().system(),
        this->dbDir(),
        *this,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    if (nccDictIO.typeHeaderOk<IOdictionary>())
    {
        IOdictionary nccDict(nccDictIO);
        if (nccDict.toc().size())
        {
            hasNonConformalCouples = true;
        }
        else if (hasNonConformalCouples)
        {
            FatalErrorInFunction
                << "Non-conformal original patches found in the mesh,"
                << " but empty 'nonConformalCouplesDict'. present"
                << exit(FatalError);
        }
    }
    else
    {
        if (hasNonConformalCouples)
        {
            FatalErrorInFunction
                << "Non-conformal original patches found in the mesh,"
                << " but cannot find file 'nonConformalCouplesDict'."
                << exit(FatalError);
        }
    }

    return hasNonConformalCouples && !stitcher_->stitches();
}


bool Foam::fvMesh::topoChanging() const
{
    return (topoChanger_.valid() ? topoChanger_->dynamic() : false);
}


bool Foam::fvMesh::dynamic() const
{
    return
    (
        (topoChanger_.valid() && topoChanger_->dynamic())
     || (mover_.valid() && mover_->dynamic())
     || (GIBChanger_.valid() && GIBChanger_->dynamic())
    );
}


bool Foam::fvMesh::update()
{
    // Note: this is a legacy function used to update the mesh for both
    // motion and topology change, which is now accessed only by a few
    // specific deprecated solvers. Please, refer to the fvMesh::topoChange
    // and fvMesh::move functions below.

    bool updated = false;

    if (GIBChanger_->dynamic())
    {
        const bool changed = GIBChanger_->update();
        topoChanged_ = changed;

        return changed;
    }

    if (!conformal()) stitcher_->disconnect(true, true);

    if (curTimeIndex_ < time().timeIndex())
    {
        const bool hasV00 = V00Ptr_;
        deleteDemandDrivenData(V00Ptr_);

        if (!hasV00)
        {
            deleteDemandDrivenData(V0Ptr_);
        }

        updated = topoChanger_->update() || updated;

        // Register V0 for distribution
        if (V0Ptr_)
        {
            V0Ptr_->checkIn();
        }

        updated = distributor_->update() || updated;

        // De-register V0 after distribution
        if (V0Ptr_)
        {
            V0Ptr_->checkOut();
        }

        if (hasV00)
        {
            // If V00 had been set reset to the mapped V0 prior to mesh-motion
            V00();
        }
    }

    updated = move() || updated;

    curTimeIndex_ = time().timeIndex();

    return updated;
}


bool Foam::fvMesh::topoChange()
{
    if (GIBChanger_->dynamic())
    {
        // Remove the oldest cell volume field
        if (V00Ptr_)
        {
            nullDemandDrivenData(V00Ptr_);
        }

        return false;
    }

    if
    (
        stitcher_->stitches()
     || topoChanger_->dynamic()
     || distributor_->dynamic()
    )
    {
        nullOldestTimeFields();
    }

    if (!conformal()) stitcher_->disconnect(true, true);

    // Remove the oldest cell volume field
    if (V00Ptr_)
    {
        nullDemandDrivenData(V00Ptr_);
    }
    else
    {
        nullDemandDrivenData(V0Ptr_);
    }

    // Set topoChanged_ false before any mesh change
    topoChanged_ = false;
    bool updated = topoChanger_->update();
    topoChanged_ = updated;

    updated = distributor_->update() || updated;

    return updated;
}


bool Foam::fvMesh::move()
{
    if (GIBChanger_->dynamic())
    {
        topoChanged_ = GIBChanger_->update();

        return topoChanged_;
    }

    if (!conformal()) stitcher_->disconnect(true, true);

    if (curTimeIndex_ < time().timeIndex() && stitcher_->stitches())
    {
        // Store all old-time fields. If we don't do this then we risk
        // triggering a store in the middle of mapping and potentially
        // overwriting a mapped old-time field with a not-yet-mapped
        // new-time field.
        storeOldTimeFields();
    }

    // Do not set moving false.
    // Once the mesh starts moving it is considered to be moving
    // for the rest of the run.
    const bool moved = mover_->update();

    curTimeIndex_ = time().timeIndex();

    stitcher_->connect(true, true, false);

    return moved;
}


void Foam::fvMesh::clearMeshPhi()
{
    // Clear mesh motion flux
    deleteDemandDrivenData(phiPtr_);
}


void Foam::fvMesh::clearV0()
{
    // Clear mesh motion flux
    deleteDemandDrivenData(V0Ptr_);
}


void Foam::fvMesh::clearGeomData()
{
    deleteDemandDrivenData(SfSlicePtr_);
    deleteDemandDrivenData(SfPtr_);
    deleteDemandDrivenData(magSfSlicePtr_);
    deleteDemandDrivenData(magSfPtr_);
    deleteDemandDrivenData(CSlicePtr_);
    deleteDemandDrivenData(CPtr_);
    deleteDemandDrivenData(CfSlicePtr_);
    deleteDemandDrivenData(CfPtr_);
}


void Foam::fvMesh::clearOut()
{
    clearGeom();

    surfaceInterpolation::clearOut();

    clearAddressing();

    polyMesh::clearOut();
}


void Foam::fvMesh::topoChange(const polyTopoChangeMap& map)
{
    // Update polyMesh. This needs to keep volume existent!
    polyMesh::topoChange(map);

    // Clear the sliced fields
    clearGeomNotOldVol();

    // Check that we're not trying to maintain old-time mesh geometry
    if (V0Ptr_ && Foam::notNull(V0Ptr_))
    {
        FatalErrorInFunction
            << "It is not possible to use mesh motion, topology change and "
            << "second-order time schemes simultaneously"
            << exit(FatalError);
    }

    // Map all fields
    mapFields(map);

    // Clear the current volume and other geometry factors
    surfaceInterpolation::clearOut();

    // Clear any non-updateable addressing
    clearAddressing(true);

    meshObject::topoChange<fvMesh>(*this, map);
    meshObject::topoChange<lduMesh>(*this, map);

    const_cast<Time&>(time()).functionObjects().topoChange(map);

    if (stitcher_.valid())
    {
        stitcher_->topoChange(map);
    }

    if (topoChanger_.valid())
    {
        topoChanger_->topoChange(map);
    }

    if (distributor_.valid())
    {
        distributor_->topoChange(map);
    }

    if (mover_.valid())
    {
        mover_->topoChange(map);
    }
}


void Foam::fvMesh::mapMesh(const polyMeshMap& map)
{
    // Distribute polyMesh data
    polyMesh::mapMesh(map);

    // Clear the sliced fields
    clearGeomNotOldVol();

    // Clear the current volume and other geometry factors
    surfaceInterpolation::clearOut();

    // Clear any non-updateable addressing
    clearAddressing(true);

    meshObject::mapMesh<fvMesh>(*this, map);
    meshObject::mapMesh<lduMesh>(*this, map);

    const_cast<Time&>(time()).functionObjects().mapMesh(map);

    stitcher_->mapMesh(map);
    topoChanger_->mapMesh(map);
    distributor_->mapMesh(map);
    mover_->mapMesh(map);
}


void Foam::fvMesh::distribute(const polyDistributionMap& map)
{
    // Distribute polyMesh data
    polyMesh::distribute(map);

    // Clear the sliced fields
    clearGeomNotOldVol();

    // Clear the current volume and other geometry factors
    surfaceInterpolation::clearOut();

    // Clear any non-updateable addressing
    clearAddressing(true);

    meshObject::distribute<fvMesh>(*this, map);
    meshObject::distribute<lduMesh>(*this, map);

    stitcher_->distribute(map);
    topoChanger_->distribute(map);
    distributor_->distribute(map);
    mover_->distribute(map);
}


void Foam::fvMesh::updateGIB()
{
    polyMesh::updateGIB();

    UpdateGIBFields<scalar, fvPatchField, volMesh>(*this);
    UpdateGIBFields<vector, fvPatchField, volMesh>(*this);
    UpdateGIBFields<sphericalTensor, fvPatchField, volMesh>(*this);
    UpdateGIBFields<symmTensor, fvPatchField, volMesh>(*this);
    UpdateGIBFields<tensor, fvPatchField, volMesh>(*this);

    UpdateGIBFields<scalar, fvsPatchField, surfaceMesh>(*this);
    UpdateGIBFields<vector, fvsPatchField, surfaceMesh>(*this);
    UpdateGIBFields<sphericalTensor, fvsPatchField, surfaceMesh>(*this);
    UpdateGIBFields<symmTensor, fvsPatchField, surfaceMesh>(*this);
    UpdateGIBFields<tensor, fvsPatchField, surfaceMesh>(*this);

    updateGeomNotOldVol();
    surfaceInterpolation::clearOut();
    clearAddressing(true);
    meshObject::movePoints<fvMesh>(*this);
    meshObject::movePoints<lduMesh>(*this);
}


void Foam::fvMesh::storeGIBFields()
{
    polyMesh::updateGIB();
    StoreGIBFields<scalar, fvPatchField, volMesh>(*this);
    StoreGIBFields<vector, fvPatchField, volMesh>(*this);
    StoreGIBFields<sphericalTensor, fvPatchField, volMesh>(*this);
    StoreGIBFields<symmTensor, fvPatchField, volMesh>(*this);
    StoreGIBFields<tensor, fvPatchField, volMesh>(*this);

    StoreGIBFields<scalar, fvsPatchField, surfaceMesh>(*this);
    StoreGIBFields<vector, fvsPatchField, surfaceMesh>(*this);
    StoreGIBFields<sphericalTensor, fvsPatchField, surfaceMesh>(*this);
    StoreGIBFields<symmTensor, fvsPatchField, surfaceMesh>(*this);
    StoreGIBFields<tensor, fvsPatchField, surfaceMesh>(*this);

    updateGeomNotOldVol();
}


void Foam::fvMesh::setPoints(const pointField& p)
{
    polyMesh::setPoints(p);

    clearGeom();

    // Update other local data
    boundary_.movePoints();
    surfaceInterpolation::movePoints();

    meshObject::movePoints<fvMesh>(*this);
    meshObject::movePoints<lduMesh>(*this);

    const_cast<Time&>(time()).functionObjects().movePoints(*this);
}


Foam::tmp<Foam::scalarField> Foam::fvMesh::movePoints(const pointField& p)
{
    // Set the mesh to be moving. This remains true for the rest of the run.
    moving_ = true;

    if (calcSolverQties())
    {
        // Create old-time volumes (if necessary) at the start of a new timestep
        if (curTimeIndex_ < time().timeIndex())
        {
            if (V00Ptr_ && notNull(V00Ptr_))
            {
                FatalErrorInFunction
                    << "Old-old volumes should not be maintained across mesh "
                    << "changes" << exit(FatalError);
            }

            // If old-old-volumes are necessary then copy them from
            // the old-volumes
            if (Foam::isNull(V00Ptr_))
            {
                V00Ptr_ = new DimensionedField<scalar, volMesh>
                (
                    IOobject
                    (
                        "Vc00",
                        time().timeName(),
                        *this,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        true
                    ),
                    V0()
                );
            }

            // Copy old-volumes from the volumes
            if (!V0Ptr_ || Foam::isNull(V0Ptr_))
            {
                V0Ptr_ = new DimensionedField<scalar, volMesh>
                (
                    IOobject
                    (
                        "Vc0",
                        time().timeName(),
                        *this,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        true
                    ),
                    V()
                );
            }
            else
            {
                V0Ptr_->scalarField::operator=(V());
            }
        }

        // Create mesh motion flux, if necessary.
        if (!phiPtr_)
        {
            phiPtr_ = new surfaceScalarField
            (
                IOobject
                (
                    "meshPhi",
                    this->time().timeName(),
                    *this,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    true
                ),
                *this,
                dimVolume/dimTime
            );
        }
        else
        {
            phiPtr_->storeOldTimes();
        }
    }

    // Move the polyMesh and set the mesh motion fluxes to the swept-volumes
    tmp<scalarField> tsweptVols = polyMesh::movePoints(p);

    if (calcSolverQties())
    {
        scalarField& sweptVols = tsweptVols.ref();
        surfaceScalarField& phi = *phiPtr_;

        scalar rDeltaT = 1.0/time().deltaTValue();

        phi.primitiveFieldRef() =
            scalarField::subField(sweptVols, nInternalFaces());
        phi.primitiveFieldRef() *= rDeltaT;

        const fvPatchList& patches = boundary();

        surfaceScalarField::Boundary& phibf = phi.boundaryFieldRef();

        forAll(patches, patchi)
        {
            phibf[patchi] = patches[patchi].patchSlice(sweptVols);
            phibf[patchi] *= rDeltaT;
        }
    }

    // Update or delete the local geometric properties as early as possible so
    // they can be used if necessary. These get recreated here instead of
    // demand driven since they might do parallel transfers which can conflict
    // with when they're actually being used.
    // Note that between above "polyMesh::movePoints(p)" and here nothing
    // should use the local geometric properties.
    updateGeomNotOldVol();

    // Update other local data
    boundary_.movePoints();
    surfaceInterpolation::movePoints();

    meshObject::movePoints<fvMesh>(*this);
    meshObject::movePoints<lduMesh>(*this);

    const_cast<Time&>(time()).functionObjects().movePoints(*this);

    return tsweptVols;
}


Foam::tmp<Foam::scalarField> Foam::fvMesh::moveGIBPoints(const pointField& p)
{
    // Set moving_ true.
    // Note: once set it remains true for the rest of the run.
    moving_ = true;

    if (V00Ptr_ && notNull(V00Ptr_))
    {
        FatalErrorInFunction
            << "Old-old volumes should not be maintained across mesh "
            << "changes" << exit(FatalError);
    }

    // If old-old-volumes are necessary then copy them from the old-volumes
    if (Foam::isNull(V00Ptr_))
    {
        V00Ptr_ = new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "Vc00",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                true
            ),
            V0()
        );
    }

    if (V0Ptr_)
    {
        // Copy V into V0 storage
        V0Ptr_->scalarField::operator=(V());
    }
    else
    {
        // Allocate V0 storage, fill with V
        V0Ptr_ = new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "Vc0",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            *this,
            dimVolume
        );
        scalarField& V0 = *V0Ptr_;
        // Note: V0 now sized with current mesh, not with (potentially
        //       different size) V.
        V0.setSize(V().size());
        V0 = V();
    }

    if (!phiPtr_)
    {
        // Create mesh motion flux
        phiPtr_ = new surfaceScalarField
        (
            IOobject
            (
                "meshPhi",
                this->time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            *this,
            dimVolume/dimTime
        );
    }
    else
    {
        // Grab old time mesh motion fluxes if the time has been incremented
        if (phiPtr_->timeIndex() != time().timeIndex())
        {
            phiPtr_->oldTime();
        }
    }

    surfaceScalarField& phi = *phiPtr_;

    // Move the polyMesh and set the mesh motion fluxes to the swept-volumes

    scalar rDeltaT = 1.0/time().deltaTValue();

    tmp<scalarField> tsweptVols = polyMesh::movePoints(p);
    scalarField& sweptVols = tsweptVols.ref();

    // Clear geometry
    clearGeomNotOldVol();

    phi.primitiveFieldRef() =
        scalarField::subField(sweptVols, nInternalFaces());
    phi.primitiveFieldRef() *= rDeltaT;

    boundary_.movePoints();

    surfaceInterpolation::clearOut();

    clearAddressing(true);

    meshObject::movePoints<fvMesh>(*this);
    meshObject::movePoints<lduMesh>(*this);
    return tsweptVols;
}


const Foam::pointField& Foam::fvMesh::oldPoints() const
{
    if (GIBChanger_.valid())
    {
        return GIBChanger_->oldPoints();
    }

    return polyMesh::oldPoints();
}


void Foam::fvMesh::conform(const surfaceScalarField& phi)
{
    // Clear the geometry fields
    clearGeomNotOldVol();

    // Clear the current volume and other geometry factors
    surfaceInterpolation::clearOut();

    // Clear any non-updateable addressing
    clearAddressing(true);

    // Modify the mesh fluxes, if necessary
    if (notNull(phi) && phiPtr_)
    {
        for (label i = 0; i <= phi.nOldTimes(); i++)
        {
            setPhi().oldTime(i).forceAssign(phi.oldTime(i));
        }
    }
}


void Foam::fvMesh::unconform
(
    const surfaceLabelField::Boundary& polyFacesBf,
    const surfaceVectorField& Sf,
    const surfaceVectorField& Cf,
    const surfaceScalarField& phi,
    const bool sync
)
{
    // Clear the geometry fields
    clearGeomNotOldVol();

    // Clear the current volume and other geometry factors
    surfaceInterpolation::clearOut();

    // Clear any non-updateable addressing
    clearAddressing(true);

    // Create non-sliced copies of geometry fields
    SfRef();
    magSfRef();
    CRef();
    CfRef();

    // Set the topology
    polyFacesBfRef().forceAssign(polyFacesBf);

    // Set the face geometry
    SfRef().forceAssign(Sf);
    magSfRef().forceAssign(mag(Sf));
    CRef().boundaryFieldRef().forceAssign(Cf.boundaryField());
    CfRef().forceAssign(Cf);

    // Communicate processor-coupled cell geometry. Cell-centre processor patch
    // fields must contain the (transformed) cell-centre locations on the other
    // side of the coupling. This is so that non-conformal patches can
    // construct weights and deltas without reference to the poly mesh
    // geometry.
    //
    // Note that the initEvaluate/evaluate communication does a transformation,
    // but it is wrong in this instance. A vector field gets transformed as if
    // it were a displacement, but the cell-centres need a positional
    // transformation. That's why there's the un-transform and re-transform bit
    // below just after the evaluate call.
    //
    // This transform handling is a bit of a hack. It would be nicer to have a
    // field attribute which identifies a field as needing a positional
    // transformation, and for it to apply automatically within the coupled
    // patch field. However, at the moment, the cell centres field is the only
    // vol-field containing an absolute position, so the hack is functionally
    // sufficient for now.
    if (sync && (Pstream::parRun() || !time().processorCase()))
    {
        volVectorField::Boundary& CBf = CRef().boundaryFieldRef();

        const label nReq = Pstream::nRequests();

        forAll(CBf, patchi)
        {
            if (isA<processorFvPatch>(CBf[patchi].patch()))
            {
                CBf[patchi].initEvaluate(Pstream::defaultCommsType);
            }
        }

        if
        (
            Pstream::parRun()
         && Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
        )
        {
            Pstream::waitRequests(nReq);
        }

        forAll(CBf, patchi)
        {
            if (isA<processorFvPatch>(CBf[patchi].patch()))
            {
                CBf[patchi].evaluate(Pstream::defaultCommsType);

                const transformer& t =
                    refCast<const processorFvPatch>(CBf[patchi].patch())
                   .transform();

                t.invTransform(CBf[patchi], CBf[patchi]);
                t.transformPosition(CBf[patchi], CBf[patchi]);
            }
        }
    }

    // Modify the mesh fluxes, if necessary
    if (notNull(phi) && phiPtr_)
    {
        for (label i = 0; i <= phi.nOldTimes(); i++)
        {
            setPhi().oldTime(i).forceAssign(phi.oldTime(i));
        }
    }

    // If debugging, do an early sanity check for the delta vectors of all
    // non-conformal cyclics. A null vector for any finite-volume face will
    // trigger an unspecified floating-point exception later on.
    if (nonConformalPolyPatch::debug)
    {
        forAll(boundary(), patchi)
        {
            const fvPatch& fvp = boundary()[patchi];

            if (!isA<nonConformalFvPatch>(fvp)) continue;

            const tmp<vectorField> tDelta(fvp.delta());
            const vectorField& delta = tDelta.ref();

            forAll(fvp, facei)
            {
                if (mag(delta[facei]) == 0)
                {
                    FatalErrorInFunction
                        << "The finite-volume face #" << facei << " of the "
                        << fvp.name() << " patch has a delta vector magnitude "
                        << "of 0," << nl << "i.e. its centre has a zero "
                        << "distance to its owner cell centroid, at point "
                        << fvp.Cf()[facei] << "." << nl << "This is an invalid "
                        << "geometry. Please check the mesh quality of the "
                        << "original patches that generated this "
                        << "non-conformal cyclic patch."
                        << exit(FatalError);
                }
            }
        }
    }
}


void Foam::fvMesh::mapFields(const polyTopoChangeMap& map)
{
    if (debug)
    {
        InfoInFunction
            << " nOldCells:" << map.nOldCells()
            << " nCells:" << nCells()
            << " nOldFaces:" << map.nOldFaces()
            << " nFaces:" << nFaces()
            << endl;
    }

    // We require geometric properties valid for the old mesh
    if
    (
        map.cellMap().size() != nCells()
     || map.faceMap().size() != nFaces()
    )
    {
        FatalErrorInFunction
            << "polyTopoChangeMap does not correspond to the old mesh."
            << " nCells:" << nCells()
            << " cellMap:" << map.cellMap().size()
            << " nOldCells:" << map.nOldCells()
            << " nFaces:" << nFaces()
            << " faceMap:" << map.faceMap().size()
            << " nOldFaces:" << map.nOldFaces()
            << exit(FatalError);
    }

    // Create a fv mapper
    const fvMeshMapper fvMap(*this, map);

    // Map all the volFields in the objectRegistry
    #define mapVolFieldType(Type, nullArg)                                     \
        MapGeometricFields<Type, fvPatchField, fvMeshMapper, volMesh>(fvMap);
    FOR_ALL_FIELD_TYPES(mapVolFieldType);

    //- Map vectorN-type volFields
    MapGeometricFields<vector1, fvPatchField, fvMeshMapper, volMesh>(fvMap);
    MapGeometricFields<vector4, fvPatchField, fvMeshMapper, volMesh>(fvMap);
    MapGeometricFields<tensor4, fvPatchField, fvMeshMapper, volMesh>(fvMap);

    // Map all the surfaceFields in the objectRegistry
    #define mapSurfaceFieldType(Type, nullArg)                                 \
        MapGeometricFields<Type, fvsPatchField, fvMeshMapper, surfaceMesh>     \
        (fvMap);
    FOR_ALL_FIELD_TYPES(mapSurfaceFieldType);

    // Map all the dimensionedFields in the objectRegistry
    #define mapVolInternalFieldType(Type, nullArg)                             \
        MapDimensionedFields<Type, fvMeshMapper, volMesh>(fvMap);
    FOR_ALL_FIELD_TYPES(mapVolInternalFieldType);

    //- Map vectorN-type dimensionedFields
    MapDimensionedFields<vector1, fvMeshMapper, volMesh>(fvMap);
    MapDimensionedFields<vector4, fvMeshMapper, volMesh>(fvMap);
    MapDimensionedFields<tensor4, fvMeshMapper, volMesh>(fvMap);

    if (pointMesh::found(*this))
    {
        // Create the pointMesh mapper
        const pointMeshMapper mapper(pointMesh::New(*this), map);

        #define mapPointFieldType(Type, nullArg)                               \
            MapGeometricFields                                                 \
            <                                                                  \
                Type,                                                          \
                pointPatchField,                                               \
                pointMeshMapper,                                               \
                pointMesh                                                      \
            >                                                                  \
            (mapper);
        FOR_ALL_FIELD_TYPES(mapPointFieldType);
    }

    // Map all the clouds in the objectRegistry
    mapClouds(*this, map);
}


void Foam::fvMesh::addPatch
(
    const label insertPatchi,
    const polyPatch& patch,
    const dictionary& patchFieldDict,
    const word& defaultPatchFieldType,
    const bool validBoundary
)
{
    // Remove my local data (see topoChange)
    // Clear mesh motion flux
    deleteDemandDrivenData(phiPtr_);

    // Clear the sliced fields
    clearGeomNotOldVol();

    // Clear the current volume and other geometry factors
    surfaceInterpolation::clearOut();

    // Clear any non-updateable addressing
    clearAddressing(true);

    const label sz = boundary_.size();

    polyMesh::addPatch
    (
        insertPatchi,
        patch,
        patchFieldDict,
        defaultPatchFieldType,
        validBoundary
    );

    boundary_.setSize(sz+1);
    boundary_.set
    (
        insertPatchi,
        fvPatch::New
        (
            boundaryMesh()[insertPatchi],
            boundary_
        )
    );

    objectRegistry& db = const_cast<objectRegistry&>(thisDb());
    AddPatchFields<volScalarField>
    (
        db,
        insertPatchi,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    AddPatchFields<volVectorField>
    (
        db,
        insertPatchi,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    AddPatchFields<volSphericalTensorField>
    (
        db,
        insertPatchi,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    AddPatchFields<volSymmTensorField>
    (
        db,
        insertPatchi,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    AddPatchFields<volTensorField>
    (
        db,
        insertPatchi,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );

    // Surface fields

    AddPatchFields<surfaceScalarField>
    (
        db,
        insertPatchi,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    AddPatchFields<surfaceVectorField>
    (
        db,
        insertPatchi,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    AddPatchFields<surfaceSphericalTensorField>
    (
        db,
        insertPatchi,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    AddPatchFields<surfaceSymmTensorField>
    (
        db,
        insertPatchi,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    AddPatchFields<surfaceTensorField>
    (
        db,
        insertPatchi,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );

    // VectorN types

    AddPatchFields<volVector1Field>
    (
        db,
        insertPatchi,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    AddPatchFields<volVector4Field>
    (
        db,
        insertPatchi,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );
    AddPatchFields<volTensor4Field>
    (
        db,
        insertPatchi,
        patchFieldDict,
        defaultPatchFieldType,
        Zero
    );

}


void Foam::fvMesh::reorderPatches
(
    const labelUList& newToOld,
    const bool validBoundary
)
{
    polyMesh::reorderPatches(newToOld, validBoundary);

    boundary_.shuffle(newToOld, validBoundary);

    objectRegistry& db = const_cast<objectRegistry&>(thisDb());

    ReorderPatchFields<volScalarField>(db, newToOld);
    ReorderPatchFields<volVectorField>(db, newToOld);
    ReorderPatchFields<volSphericalTensorField>(db, newToOld);
    ReorderPatchFields<volSymmTensorField>(db, newToOld);
    ReorderPatchFields<volTensorField>(db, newToOld);

    ReorderPatchFields<surfaceScalarField>(db, newToOld);
    ReorderPatchFields<surfaceVectorField>(db, newToOld);
    ReorderPatchFields<surfaceSphericalTensorField>(db, newToOld);
    ReorderPatchFields<surfaceSymmTensorField>(db, newToOld);
    ReorderPatchFields<surfaceTensorField>(db, newToOld);

    ReorderPatchFields<volVector1Field>(db, newToOld);
    ReorderPatchFields<volVector4Field>(db, newToOld);
    ReorderPatchFields<volTensor4Field>(db, newToOld);
}


void Foam::fvMesh::removeFvBoundary()
{
    if (debug)
    {
        InfoInFunction << "Removing boundary patches." << endl;
    }

    // Remove fvBoundaryMesh data first.
    boundary_.clear();
    boundary_.setSize(0);
    polyMesh::removeBoundary();

    clearOut();
}


void Foam::fvMesh::swap(fvMesh& otherMesh)
{
    // Clear the sliced fields
    clearGeom();

    // Clear the current volume and other geometry factors
    surfaceInterpolation::clearOut();

    // Clear any non-updateable addressing
    clearAddressing(true);

    polyMesh::swap(otherMesh);

    auto updatePatches = []
    (
        const polyPatchList& patches,
        fvBoundaryMesh& boundaryMesh
    )
    {
        boundaryMesh.setSize(patches.size());

        forAll(patches, patchi)
        {
            // Construct new processor patches, as the decomposition may have
            // changed. Leave other patches as is.

            if (isA<processorPolyPatch>(patches[patchi]))
            {
                boundaryMesh.set
                (
                    patchi,
                    fvPatch::New
                    (
                        patches[patchi],
                        boundaryMesh
                    )
                );
            }
        }
    };

    updatePatches(boundaryMesh(), boundary_);
    updatePatches(otherMesh.boundaryMesh(), otherMesh.boundary_);
}


Foam::tmp<Foam::vectorField> Foam::fvMesh::velocityCorrect
(
    const vectorField& pc
) const
{
    tmp<vectorField> zeroFieldt
        (
            new vectorField
            (
                pc.size(), vector::zero
            )
        );
    return zeroFieldt;
}


bool Foam::fvMesh::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    bool ok = true;

    if (!conformal() && pointsWriteOpt() == IOobject::AUTO_WRITE)
    {
        // Create a full surface field with the polyFacesBf boundary field to
        // write to disk. Make the internal field uniform to save disk space.

        surfaceLabelField polyFaces
        (
            polyFacesBfIO(IOobject::NO_READ),
            *this,
            dimless,
            labelField(nInternalFaces(), -1),
            *polyFacesBfPtr_
        );

        ok = ok & polyFaces.write(write);
    }

    if (phiPtr_)
    {
        ok = ok && phiPtr_->write(write);

        // NOTE: The old old time mesh phi might be necessary for certain
        // solver smooth restart using second order time schemes.
        //ok = phiPtr_->oldTime().write(write);
    }

    // For second-order restarts we need to write V0
    if (V00Ptr_)
    {
        ok = ok && V0Ptr_->write(write);
    }

    if (stitcher_.valid())
    {
        stitcher_->write(write);
    }

    if (topoChanger_.valid())
    {
        topoChanger_->write(write);
    }

    if (distributor_.valid())
    {
        distributor_->write(write);
    }

    if (mover_.valid())
    {
        mover_->write(write);
    }

    if (GIBChanger_.valid())
    {
        GIBChanger_->write(write);
    }

    return ok && polyMesh::writeObject(fmt, ver, cmp, write);
}


bool Foam::fvMesh::write(const bool write) const
{
    return polyMesh::write(write);
}


const Foam::fvSchemes& Foam::fvMesh::schemes() const
{
    if (!fvSchemes_.valid())
    {
        fvSchemes_ = new fvSchemes(*this);
    }

    return fvSchemes_();
}


const Foam::fvSolution& Foam::fvMesh::solution() const
{
    if (!fvSolution_.valid())
    {
        fvSolution_ = new fvSolution(*this);
    }

    return fvSolution_();
}


Foam::fvSchemes& Foam::fvMesh::schemes()
{
    if (!fvSchemes_.valid())
    {
        fvSchemes_ = new fvSchemes(*this);
    }

    return *fvSchemes_;
}


Foam::fvSolution& Foam::fvMesh::solution()
{
    if (!fvSolution_.valid())
    {
        fvSolution_ = new fvSolution(*this);
    }

    return *fvSolution_;
}


template<>
//typename Foam::pTraits<Foam::sphericalTensor>::labelType
Foam::pTraits<Foam::sphericalTensor>::labelType
Foam::fvMesh::validComponents<Foam::sphericalTensor>() const
{
    return Foam::pTraits<Foam::sphericalTensor>::labelType(1);
}


template<>
Foam::pTraits<Foam::sphericalTensor>::labelType
Foam::fvMesh::validComponents2<Foam::sphericalTensor>() const
{
    return Foam::pTraits<Foam::sphericalTensor>::labelType(1);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::fvMesh::operator!=(const fvMesh& bm) const
{
    return &bm != this;
}


bool Foam::fvMesh::operator==(const fvMesh& bm) const
{
    return &bm == this;
}


// ************************************************************************* //
