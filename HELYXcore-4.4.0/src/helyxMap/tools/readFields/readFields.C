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
    (c) 2019-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "tools/readFields/readFields.H"
#include "helyxMap.H"
#include "db/IOobjectList/IOobjectList.H"
#include "fields/fvsPatchFields/fvsPatchField/fvsPatchFields.H"
#include "meshes/polyMesh/polyPatches/constraint/symmetryPlane/symmetryPlanePolyPatch.H"
#include "nonConformal/polyPatches/nonConformal/nonConformalPolyPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::readFields::excludePatch(const label patchi) const
{
    const polyPatch& pp = mesh_->boundaryMesh()[patchi];

    if
    (
        isA<emptyPolyPatch>(pp)
     || isA<symmetryPlanePolyPatch>(pp)
     || isA<nonConformalPolyPatch>(pp)
    )
    {
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::readFields::readFields()
:
    type_("source"),
    norm_(false)
{}


Foam::readFields::readFields
(
    const fvMesh* mesh,
    const Time* runTime,
    helyxMap* map,
    const word& type
)
:
    type_(type),
    mapTime_(runTime->timeName()),
    mesh_(mesh),
    runTime_(runTime),
    fieldMap_(map),
    norm_(false),
    alphaMax_(fieldMap_->alphaMax_),
    nCells_(mesh->C().size()),
    scFields_(fieldMap_->mapScalarFields_),
    vecFields_(fieldMap_->mapVectorFields_),
    rhoRef_(1),
    UrotDegreeFromSource_ (0),
    Uref_(point(1, 0, 0)),
    interp_(fieldMap_->interp_)
{
    if (fieldMap_->mapTimeName_.length() >= 1)
    {
        mapTime_ = fieldMap_->mapTimeName_;
    }

    boundBox meshBb(mesh->points(), true);
    getBoundBox(meshBb);
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::readFields::~readFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::readFields::createFields(const word& tName)
{
    word timeName(tName);
    if (timeName.length() < 1)
    {
        timeName = runTime_->timeName();
    }

    IOobjectList objects(*mesh_, timeName);

    HashSet<word> removeScalars;
    HashSet<word> removeVectors;

    forAllConstIter(HashSet<word>, fieldMap_->mapScalarFields_, iter)
    {
        const word& scName = iter();

        Info<< "Reading scalar field: " << scName;

        if (objects.lookup(scName) == nullptr)
        {
            Info<< " ... not found!" << endl;

            removeScalars.insert(scName);

            continue;
        }

        Info<< endl;

        autoPtr<volScalarField> scFieldPtr
        (
            new volScalarField
            (
                IOobject
                (
                    scName,
                    timeName,
                    *mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                *mesh_
            )
        );

        scFieldPtr.ptr()->store();
    }

    forAllConstIter(HashSet<word>, fieldMap_->mapVectorFields_, iter)
    {
        const word& vecName = iter();

        Info<< "Reading vector field: " << vecName;

        if (objects.lookup(vecName) == nullptr)
        {
            Info<< " ... not found!" << endl;

            removeVectors.insert(vecName);

            continue;
        }

        Info<< endl;

        autoPtr<volVectorField> vecFieldPtr
        (
            new volVectorField
            (
                IOobject
                (
                    vecName,
                    timeName,
                    *mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                *mesh_
            )
        );

        vecFieldPtr.ptr()->store();
    }

    fieldMap_->mapScalarFields_ -= removeScalars;
    fieldMap_->mapVectorFields_ -= removeVectors;

    // Flux map
    forAllConstIter(HashSet<word>, fieldMap_->mapSurfaceScalarFields_, iter)
    {
        const word& scName = iter();

        autoPtr<surfaceScalarField> scFieldPtr
        (
            new surfaceScalarField
            (
                IOobject
                (
                    scName,
                    timeName,
                    *mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                *mesh_
            )
        );

        scFieldPtr.ptr()->store();
    }
}


void Foam::readFields::storeFields()
{
    if (type_ != "source") return;

    label nBFace = numBoundFaces();

    label nCell = mesh_->C().size();
    label sz = nCell;

    if (fieldMap_->mapBoundary_)
    {
        sz += nBFace;
    }

    xyz_.setSize(sz);
    forAll(mesh_->C(), i)
    {
        xyz_[i] = mesh_->C()[i];
    }

    if (fieldMap_->mapBoundary_)
    {
        label cnt = nCell;
        forAll(mesh_->boundaryMesh(), r)
        {
            if (excludePatch(r)) continue;

            const vectorField& bndPoints = mesh_->Cf().boundaryField()[r];
            forAll(bndPoints, j)
            {
                const point& pt = bndPoints[j];
                xyz_[cnt] = pt;
                cnt++;
            }
        }
    }

    forAllConstIter(HashSet<word>, scFields_, iter)
    {
        const word& name = iter();
        scalar coef = scale(name);

        const volScalarField& sField = getScalarField(name);
        label sz = sField.primitiveField().size();
        if (fieldMap_->mapBoundary_)
        {
            sz += nBFace;
        }

        scalarFields_[name].resize(sz);

        for (label i = 0; i < nCell; i++)
        {
            scalarFields_[name][i] = sField.primitiveField()[i]/coef;
        }

        if (fieldMap_->mapBoundary_)
        {
            // Boundaries
            label cnt = nCell;
            forAll(mesh_->boundaryMesh(), r)
            {
                if (excludePatch(r)) continue;

                forAll(mesh_->boundaryMesh()[r], j)
                {
                    scalar sc = sField.boundaryField()[r][j];
                    scalarFields_[name][cnt] = sc/coef;
                    cnt++;
                }
            }
        }
    }

    // Vector fields
    forAllConstIter(HashSet<word>, vecFields_, iter)
    {
        const word& name = iter();
        scalar coef = scale(name);

        const volVectorField& vsField = getVectorField(name);

        label sz = vsField.primitiveField().size();
        if (fieldMap_->mapBoundary_)
        {
            sz += nBFace;
        }

        vectorFields_[name].resize(sz);
        for (label i = 0; i < nCell; i++)
        {
            vectorFields_[name][i] = vsField.primitiveField()[i]/coef;
        }

        if (fieldMap_->mapBoundary_)
        {
            // Boundaries
            label cnt = nCell;
            forAll(mesh_->boundaryMesh(), r)
            {
                if (excludePatch(r)) continue;

                forAll(mesh_->boundaryMesh()[r], j)
                {
                    point vsc = vsField.boundaryField()[r][j];
                    vectorFields_[name][cnt] = vsc/coef;
                    cnt++;
                }
            }
        }
    }
}


void Foam::readFields::storeFields(const boundBox& sampleBox)
{
    if (type_ != "source") return;

    xyz_.resize(0);

    std::vector<label> inCells;

    typedef Tuple2<label, label> regFace;
    std::vector<regFace> bFaces;

    forAll(mesh_->C(), i)
    {
        const point& pt = mesh_->C()[i];

        if (sampleBox.contains(pt))
        {
            xyz_.append(pt);
            inCells.push_back(i);
        }
    }

    if (fieldMap_->mapBoundary_)
    {
        forAll(mesh_->boundaryMesh(), r)
        {
            if (excludePatch(r)) continue;

            const vectorField& bndPoints = mesh_->Cf().boundaryField()[r];
            forAll(bndPoints, j)
            {
                const point& pt = bndPoints[j];

                if (sampleBox.contains(pt))
                {
                    xyz_.append(pt);
                    regFace face(r, j);
                    bFaces.push_back(face);
                }
            }
        }
    }

    // Scalar fields
    forAllConstIter(HashSet<word>, scFields_, iter)
    {
        const word& name = iter();
        scalar coef = scale(name);

        // Internal field
        const volScalarField& sField = getScalarField(name);

        label nCell = inCells.size();
        label nBFace = bFaces.size();

        scalarFields_[name].resize(nCell + nBFace);

        for (label i = 0; i < nCell; i++)
        {
            label ic = inCells[i];
            scalarFields_[name][i] = sField.primitiveField()[ic]/coef;
        }

        // Boundaries
        if (fieldMap_->mapBoundary_)
        {
            for (label jf = 0; jf < nBFace; jf++)
            {
                const regFace& face = bFaces[jf];

                label r = face.first();
                label j = face.second();

                scalar sc = sField.boundaryField()[r][j];

                scalarFields_[name][nCell + jf] = sc/coef;
            }
        }
    }

    // Vector fields
    forAllConstIter(HashSet<word>, vecFields_, iter)
    {
        const word& name = iter();
        scalar coef = scale(name);

        // Internal field
        const volVectorField& vsField = getVectorField(name);

        label nCell = inCells.size();
        label nBFace = bFaces.size();

        vectorFields_[name].resize(nCell + nBFace);

        for (label i = 0; i < nCell; i++)
        {
            label ic = inCells[i];

            vectorFields_[name][i] = vsField.primitiveField()[ic]/coef;
        }

        // Boundaries
        if (fieldMap_->mapBoundary_)
        {
            for (label jf = 0; jf < nBFace; jf++)
            {
                const regFace& face = bFaces[jf];

                label r = face.first();
                label j = face.second();

                point vsc = vsField.boundaryField()[r][j];

                vectorFields_[name][nCell + jf] = vsc/coef;
            }
        }
    }
}


void Foam::readFields::setWallDist(const volScalarField& y, scalar maxWDist)
{
    maxWallDist_ = maxWDist;

    label nCells = y.primitiveField().size();
    label nBFace = numBoundFaces();

    label sz = nCells;
    if (fieldMap_->mapBoundary_)
    {
        sz += nBFace;
    }

    wDists_.setSize(sz);
    wDists_ = 1;

    for (label i = 0; i < nCells; i++)
    {
        wDists_[i] = y.primitiveField()[i]/maxWallDist_;
    }

    deltaY_ = 1.0/float(fieldMap_->nwdist_);
}


void Foam::readFields::buildKdTrees
(
    const volScalarField& y,
    scalar cutoff,
    const word& excludeBndType
)
{
    maxWallDist_ = gMax(y.primitiveField());

    label nCells = y.primitiveField().size();

    wDists_.setSize(nCells);
    for (label i = 0; i < nCells; i++)
    {
        wDists_[i] = y.primitiveField()[i]/maxWallDist_;
    }

    deltaY_ = 1.0/float(fieldMap_->nwdist_);

    Array2d<point> xyz;
    xyz.resize(fieldMap_->nwdist_ + 1);

    alphaMap_.resize(fieldMap_->nwdist_ + 1);

    forAll(mesh_->C(), i)
    {
        point pt = mesh_->C()[i];

        if (wDists_[i] > fieldMap_->alphaMax_) continue;

        label ic = inGroup(i);
        if (ic < 0)
        {
            ic = 0;
        }
        if (ic > fieldMap_->nwdist_)
        {
            ic = fieldMap_->nwdist_;
        }

        scalar xb = xbar(pt.x());
        scalar yb = ybar(pt.y());
        scalar zb = zbar(pt.z());

        point pt1(xb, yb, zb);

        xyz[ic].push_back(pt1);
        alphaMap_[ic].push_back(i);
    }

    if (type_ == "source")
    {
        kdTrees_.resize(fieldMap_->nwdist_ + 1);

        for (unsigned ic = 0; ic < xyz.size(); ic++)
        {
            label npt = xyz[ic].size();
            if (npt < 1) continue;

            farray2d trainData(npt);
            for (label i = 0; i < npt; i++)
            {
                trainData[i].resize(3);

                point pt = xyz[ic][i];
                trainData[i][0] = pt.x();
                trainData[i][1] = pt.y();
                trainData[i][2] = pt.z();
            }

            kdTrees_[ic].build(trainData);
        }

        // Boundary search tree
        constructInternalKnn(cutoff);

        if (fieldMap_->mapBoundary_)
        {
            constructBoundaryKnn();
        }
    }
}


void Foam::readFields::getBoundBox(const boundBox& box)
{
    const point& pmin = box.min();
    const point& pmax = box.max();

    xmin_ = pmin.x();
    ymin_ = pmin.y();
    zmin_ = pmin.z();
    xmax_ = pmax.x();
    ymax_ = pmax.y();
    zmax_ = pmax.z();
}


void Foam::readFields::setInputs()
{
    if (type_ == "source")
    {
        Uref_ = fieldMap_->UrefSource_;
        rhoRef_ = fieldMap_->rhoRefSource_;
    }
    else
    {
        Uref_ = fieldMap_->UrefTarget_;
        rhoRef_ = fieldMap_->rhoRefTarget_;
    }

    interp_ = fieldMap_->interp_;
}


Foam::scalar Foam::readFields::scale(const word& fldName)
{
    const word unknown("unknown");
    const word& fType = fieldMap_->fieldTypes_.lookup(fldName, unknown);

    if (fType == "velocity")
    {
        return U0();
    }
    else if (fType == "pressure")
    {
        return pRef();
    }
    else if (fType == "turbEnergy")
    {
        return U0()*U0();
    }

    return 1;
}


Foam::label Foam::readFields::numBoundFaces()
{
    label nBFace = 0;

    forAll(mesh_->boundaryMesh(), r)
    {
        if (excludePatch(r)) continue;

        label sz = mesh_->boundaryMesh()[r].size();
        nBFace += sz;
    }

    return nBFace;
}


void Foam::readFields::constructInternalKnn(scalar cutoff)
{
    Info<< "\nConstructing internal KNN map..." << endl;

    label cnt = 0;
    cutoff = -0.1;

    farray2d trainData;

    // For gridMap (cutoff < 0), no need to do alpha-mapping,
    // internal map contains all cells.
    forAll(mesh_->C(), i)
    {
        if (cutoff >= 0)
        {
            scalar alphai = wDists_[i];

            if (alphai > cutoff) continue;
        }

        point pt = mesh_->C()[i];
        scalar xb = xbar(pt.x());
        scalar yb = ybar(pt.y());
        scalar zb = zbar(pt.z());

        std::vector<scalar> pts(3);
        pts[0] = xb;
        pts[1] = yb;
        pts[2] = zb;

        trainData.push_back(pts);

        cellMap_[cnt] = i;
        cnt++;
    }

    Info<< "Train data size: " << trainData.size() << " " << cnt << endl;
    internalMap_.build(trainData);

    Info<< "Internal map build." << endl;
}


void Foam::readFields::constructBoundaryKnn()
{
    farray2d trainData;

    forAll(mesh_->boundaryMesh(), r)
    {
        if (excludePatch(r)) continue;

        word name = mesh_->boundaryMesh()[r].name();

        forAll(mesh_->boundaryMesh()[r], j)
        {
            point pt = mesh_->Cf().boundaryField()[r][j];

            scalar xb = xbar(pt.x());
            scalar yb = ybar(pt.y());
            scalar zb = zbar(pt.z());

            std::vector<scalar> pts(3);
            pts[0] = xb;
            pts[1] = yb;
            pts[2] = zb;

            trainData.push_back(pts);
        }
    }

    if (trainData.size() >= 1)
    {
        boundaryMap_.build(trainData);
    }
}


void Foam::readFields::transformPoints
(
    const boundBox& tBox,
    DynamicList<point>& xyz
)
{
    // Called after construct boundBox, note that tBox must be the original
    // (not extended).
    if (myType() != "source")
    {
        return;
    }

    scalar xmint = tBox.min().x();
    scalar xmaxt = tBox.max().x();
    scalar xmins = gBox().min().x();
    scalar xmaxs = gBox().max().x();
    scalar xcoef = (xmaxt - xmint)/(xmaxs - xmins);

    scalar ymint = tBox.min().y();
    scalar ymaxt = tBox.max().y();
    scalar ymaxs = gBox().max().y();
    scalar ymins = gBox().min().y();
    scalar ycoef = (ymaxt - ymint)/(ymaxs - ymins);

    scalar zmint = tBox.min().z();
    scalar zmaxt = tBox.max().z();
    scalar zmaxs = gBox().max().z();
    scalar zmins = gBox().min().z();
    scalar zcoef = (zmaxt - zmint)/(zmaxs - zmins);

    forAll(xyz, i)
    {
        scalar x = xyz[i].x();
        x = xmint + (x - xmins)*xcoef;
        xyz[i].x() = x;

        scalar y = xyz[i].y();
        y = ymint + (y - ymins)*ycoef;
        xyz[i].y() = y;

        scalar z = xyz[i].z();
        z = zmint + (z - zmins)*zcoef;
        xyz[i].z() = z;
    }
}


void Foam::readFields::constructBoundBox
(
    scalar extCoef,
    const word type,
    const label numCells
)
{
    // Global box
    gBox_ = new boundBox(mesh_->C(), true);

    // Local box
    boundBox lBox(mesh_->C(), false);

    if (extCoef > SMALL)
    {
        lBox.inflate(extCoef);
    }

    bBox_ = new smartBoundBox(lBox.min(), lBox.max(), type);

    if (type == "smartBoundBox")
    {
        bBox_().createBox(mesh_, numCells);
    }
}


// ************************************************************************* //
