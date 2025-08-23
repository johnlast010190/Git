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

#include "surfaceMap/surfaceMap.H"
#include "fields/fvsPatchFields/fvsPatchField/fvsPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceMap::surfaceMap()
:
    helyxMap(),
    imap_(true),
    bmap_(true),
    byType_(false),
    byPhysicalType_(false),
    byNameList_(true),
    bodySurfaceType_("wall")
{}


Foam::surfaceMap::surfaceMap(const fvMesh* mesh, const Time* runTime)
:
    helyxMap(),
    mesh_(mesh),
    runTime_(runTime),
    imap_(true),
    bmap_(true),
    byType_(false),
    byPhysicalType_(false),
    byNameList_(true),
    bodySurfaceType_("wall")
{}


Foam::surfaceMap::surfaceMap(const fvMesh& mesh)
:
    helyxMap(),
    mesh_(&mesh),
    runTime_(&mesh.time()),
    imap_(true),
    bmap_(true),
    byType_(false),
    byPhysicalType_(false),
    byNameList_(true),
    bodySurfaceType_("wall")
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::surfaceMap::~surfaceMap()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfaceMap::setInput(const dictionary& dict)
{
    helyxMap::setInput(dict);

    word identifyObject =
        dict.lookupOrDefault<word>("identifyObject", "byNameList");

    if (identifyObject == "byType")
    {
        byType_ = true;
    }
    else if (identifyObject == "byPhysicalType")
    {
        byPhysicalType_ = true;
    }
    else
    {
        byNameList_ = true;
    }

    bodySurfaceType_ = dict.lookupOrDefault<word>("bodySurfaceType", "wall");

    imap_ = dict.lookupOrDefault<bool>("internalSym", true);
    bmap_ = dict.lookupOrDefault<bool>("surfaceSym", false);

    point p0(0, 0, 0);
    vector norm(0, 1, 0);

    plane plane0(p0, norm);
    plane newPlane = dict.lookupOrDefault<plane>("mirror", plane0);

    mirror_ = new plane(newPlane);

    if (byType_)
    {
        sourcePatchNames_.resize(0);
        forAll(mesh_->boundaryMesh(), r)
        {
            word type = mesh_->boundaryMesh()[r].type();
            word name = mesh_->boundaryMesh()[r].name();

            if (type == bodySurfaceType_)
            {
                sourcePatchNames_.push_back(name);
            }
        }
    }
    else if (byPhysicalType_)
    {
        sourcePatchNames_.resize(0);
        forAll(mesh_->boundaryMesh(), r)
        {
            word type = mesh_->boundaryMesh()[r].physicalType();
            word name = mesh_->boundaryMesh()[r].name();

            if (type == bodySurfaceType_)
            {
                sourcePatchNames_.push_back(name);
            }
        }
    }
    else
    {
        const polyBoundaryMesh& bMesh = mesh_->boundaryMesh();

        labelHashSet includePatches;

        if (dict.found("patches"))
        {
            includePatches = bMesh.patchSet(wordReList(dict.lookup("patches")));
        }

        const wordList allPatchNames(bMesh.names());

        forAllConstIter(labelHashSet, includePatches, it)
        {
            sourcePatchNames_.push_back(allPatchNames[*it]);
        }
    }

    Info<< "Number of selected patch regions: " << sourcePatchNames_.size()
        << endl;
}


void Foam::surfaceMap::setInput(const dictionary& dict, const argList& args)
{
    helyxMap::setInput(dict, args);

    word identifyObject =
        dict.lookupOrDefault<word>("identifyObject", "byNameList");

    if (identifyObject == "byType")
    {
        byType_ = true;
    }
    else if (identifyObject == "byPhysicalType")
    {
        byPhysicalType_ = true;
    }
    else
    {
        byNameList_ = true;
    }

    bodySurfaceType_ = dict.lookupOrDefault<word>("bodySurfaceType", "wall");

    imap_ = dict.lookupOrDefault<bool>("internalSym", true);
    bmap_ = dict.lookupOrDefault<bool>("surfaceSym", false);

    point p0(0, 0, 0);
    vector norm(0, 1, 0);

    plane plane0(p0, norm);
    plane newPlane = dict.lookupOrDefault<plane>("mirror", plane0);

    mirror_ = new plane(newPlane);

    if (byType_)
    {
        sourcePatchNames_.resize(0);
        forAll(mesh_->boundaryMesh(), r)
        {
            word type = mesh_->boundaryMesh()[r].type();
            word name = mesh_->boundaryMesh()[r].name();

            if (type == bodySurfaceType_)
            {
                sourcePatchNames_.push_back(name);
            }
        }
    }
    else if (byPhysicalType_)
    {
        sourcePatchNames_.resize(0);
        forAll(mesh_->boundaryMesh(), r)
        {
            word type = mesh_->boundaryMesh()[r].physicalType();
            word name = mesh_->boundaryMesh()[r].name();

            if (type == bodySurfaceType_)
            {
                sourcePatchNames_.push_back(name);
            }
        }
    }
    else
    {
        const polyBoundaryMesh& bMesh = mesh_->boundaryMesh();

        labelHashSet includePatches;

        if (dict.found("patches"))
        {
            includePatches = bMesh.patchSet(wordReList(dict.lookup("patches")));
        }

        const wordList allPatchNames(bMesh.names());

        forAllConstIter(labelHashSet, includePatches, it)
        {
            sourcePatchNames_.push_back(allPatchNames[*it]);
        }
    }

    Info<< "Number of selected patch regions: " << sourcePatchNames_.size()
        << endl;
}


void Foam::surfaceMap::getParallelFields()
{
    if (!Pstream::parRun())
    {
        return;
    }

    List<vectorField> posLoc(Pstream::nProcs());
    label myid = Pstream::myProcNo();

    if (internalSym())
    {
        DynamicList<point> xyz;
        xyz.resize(mesh_->C().size());
        forAll(mesh_->C(), i)
        {
            xyz[i] = mesh_->C()[i];
        }

        posLoc[myid] = xyz;

        Pstream::allGatherList(posLoc);
        vectorField posGlob =
            ListListOps::combine<vectorField>(posLoc, accessOp<vectorField>());

        sourceXyz_ = posGlob;
        posGlob.resize(0);
    }

    Info<< "Source_xyz: " << sourceXyz_.size() << endl;

    if (bndSym())
    {
        // Boundary
        List<vectorField> posLoc(Pstream::nProcs());

        for (unsigned ir = 0; ir < sourcePatchNames_.size(); ir++)
        {
            const word& bname = sourcePatchNames_[ir];
            label r = patchId(bname);
            if (r < 0)
            {
                 WarningInFunction
                    << "Warning: patch ID " << bname << " not found.\n";

                 continue;
            }

            DynamicList<point> bxyz;
            forAll(mesh_->Cf().boundaryField()[r], j)
            {
                const point& pt = mesh_->Cf().boundaryField()[r][j];
                bxyz.append(pt);
            }

            posLoc[myid] = bxyz;
            Pstream::allGatherList(posLoc);

            vectorField posGlob =
                ListListOps::combine<vectorField>
                (
                    posLoc,
                    accessOp<vectorField>()
                );

            bSourceXYZ_[bname] = posGlob;
        }
    }

    forAllConstIter(HashSet<word>, mapScalarFields_, iter)
    {
        const word& scname = iter();

        const volScalarField& fld = source().getScalarField(scname);
        label sz = fld.primitiveField().size();

        if (internalSym())
        {
            scalarField scfldLoc(sz);
            forAll(fld.primitiveField(), i)
            {
                scfldLoc[i] = fld[i];
            }

            List<scalarField> sField(Pstream::nProcs());
            sField[myid] = scfldLoc;
            Pstream::allGatherList(sField);

            scalarField fldGlobal =
                ListListOps::combine<scalarField>
                (
                    sField,
                    accessOp<scalarField>()
                );

            sourceScalarFields_[scname] = fldGlobal;

            wait();
        }

        if (bndSym())
        {
            for (unsigned ir = 0; ir < sourcePatchNames_.size(); ir++)
            {
                word pname = sourcePatchNames_[ir];
                word key = pname + ',' + scname;

                // Local field on boundary
                List<scalarField> bSField(Pstream::nProcs());

                label pid = patchId(pname);

                DynamicList<scalar> sclist;
                forAll(fld.boundaryField()[pid], j)
                {
                    scalar scb = fld.boundaryField()[pid][j];
                    sclist.append(scb);
                }

                bSField[myid] = sclist;
                Pstream::allGatherList(bSField);

                scalarField bfldGlobal =
                    ListListOps::combine<scalarField>
                    (
                        bSField,
                        accessOp<scalarField>()
                    );

                bSourceScalarFields_[key] = bfldGlobal;

                bfldGlobal.resize(0);
            }
        }
    }
}


void Foam::surfaceMap::createMirrorFields(const word& timeName)
{
    forAllConstIter(HashSet<word>, mapScalarFields_, iter)
    {
        const word& scname = iter();
        word mScalar = mirrorName(scname);
        volScalarField& fld = source().getScalarField(scname);

        autoPtr<volScalarField> scfieldPtr
        (
            new volScalarField
            (
                IOobject
                (
                    mScalar,
                    timeName,
                    *mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                *mesh_,
                fld.dimensions()
            )
        );

        scfieldPtr.ptr()->store();

        volScalarField& mfld = source().getScalarField(mScalar);
        mfld = fld;
    }
}


void Foam::surfaceMap::buildSearchTreePar()
{
    if (internalSym())
    {
        farray2d trainData;

        forAll(sourceXyz_, i)
        {
            const point& pt = sourceXyz_[i];

            scalar xb = source().xbar(pt.x());
            scalar yb = source().ybar(pt.y());
            scalar zb = source().zbar(pt.z());

            std::vector<scalar> pts(3);
            pts[0] = xb;
            pts[1] = yb;
            pts[2] = zb;

            trainData.push_back(pts);
        }

        internalMap_.build(trainData);

        wait();
    }

    if (bndSym())
    {
        std::map<word, DynamicList<point>>::iterator it;

        for (it = bSourceXYZ_.begin(); it != bSourceXYZ_.end(); it++)
        {
            word bName = it->first;
            const DynamicList<point>& bfList = it->second;

            farray2d tData;

            for (label j = 0; j < bfList.size(); j++)
            {
                const point& pt = bfList[j];

                std::vector<scalar> pts(3);
                pts[0] = pt.x();
                pts[1] = pt.y();
                pts[2] = pt.z();

                tData.push_back(pts);
                faceIdMaps_[bName].push_back(j);
            }

            bndMaps_[bName].build(tData);
        }

        wait();
    }
}


void Foam::surfaceMap::buildSearchTree()
{
    if (Pstream::parRun())
    {
        buildSearchTreePar();

        return;
    }

    if (internalSym())
    {
        farray2d trainData;

        forAll(mesh_->C(), i)
        {
            const point& pt = mesh_->C()[i];

            scalar xb = source().xbar(pt.x());
            scalar yb = source().ybar(pt.y());
            scalar zb = source().zbar(pt.z());

            std::vector<scalar> pts(3);
            pts[0] = xb;
            pts[1] = yb;
            pts[2] = zb;

            trainData.push_back(pts);
        }

        Info<< "Number of data samples: " << trainData.size() << endl;
        internalMap_.build(trainData);
    }

    if (bndSym())
    {
        for (unsigned ir = 0; ir < sourcePatchNames_.size(); ir++)
        {
            const word& bname = sourcePatchNames_[ir];

            label r = patchId(bname);
            if (r < 0)
            {
                 WarningInFunction
                    << "Warning, patch ID: "<< bname << " not found.\n";

                 continue;
            }

            farray2d tData;
            std::vector<label> faceIds;

            forAll(mesh_->Cf().boundaryField()[r], j)
            {
                const point& pt = mesh_->Cf().boundaryField()[r][j];

                std::vector<scalar> pts(3);
                pts[0] = pt.x();
                pts[1] = pt.y();
                pts[2] = pt.z();

                tData.push_back(pts);
                faceIdMaps_[bname].push_back(j);
            }

            bndMaps_[bname].build(tData);
        }

        Info<< "Number of boundary maps: " << bndMaps_.size() << endl;
    }
}


Foam::label Foam::surfaceMap::patchId(const word& name)
{
    forAll(mesh_->boundaryMesh(), r)
    {
        if (mesh_->boundaryMesh()[r].name() == name)
        {
            return r;
        }
    }

    return -1;
}


void Foam::surfaceMap::getSymFields()
{
    Info<< "Obtain the symmetric fields..." << endl;

    forAllConstIter(HashSet<word>, mapScalarFields_, iter)
    {
        const word& scName = iter();

        const volScalarField& scField = source().getScalarField(scName);

        word mscName = mirrorName(scName);

        volScalarField& mscField = source().getScalarField(mscName);
        mscField = 0.5*(mscField + scField);
    }
}


void Foam::surfaceMap::getMirrorFields()
{
    std::vector<label> uList;
    std::map<word, std::vector<label>> uBList;

    if (internalSym())
    {
        uList.resize(mesh_->C().size());

        forAll(mesh_->C(), i)
        {
            const point& pt = mesh_->C()[i];
            point mpt = mirror().mirror(pt);

            scalar xc = source().xbar(mpt.x());
            scalar yc = source().ybar(mpt.y());
            scalar zc = source().zbar(mpt.z());

            std::vector<scalar> qvect(3);
            qvect[0] = xc;
            qvect[1] = yc;
            qvect[2] = zc;

            std::vector<label> nb;
            std::vector<scalar> dist;

            internalMap_.search(qvect, nb, dist);

            label ic = nb[0];
            uList[i] = ic;
        }
    }

    if (bndSym())
    {
        for (unsigned r = 0; r < sourcePatchNames_.size(); r++)
        {
            const word& pname = sourcePatchNames_[r];

            label ir = patchId(pname);
            if (ir < 0) continue;

            forAll(mesh_->boundaryMesh()[ir], j)
            {
                const point& pt = mesh_->Cf().boundaryField()[ir][j];
                point mpt = mirror().mirror(pt);

                std::vector<scalar> pts(3);
                pts[0] = mpt.x();
                pts[1] = mpt.y();
                pts[2] = mpt.z();

                std::vector<label> nb;
                std::vector<scalar> dist;

                bndMaps_[pname].search(pts, nb, dist);

                label ib = nb[0];
                label jc = faceIdMaps_[pname][ib];

                uBList[pname].push_back(jc);
            }
        }
    }

    wait();

    forAllConstIter(HashSet<word>, mapScalarFields_, iter)
    {
        const word& scName = iter();
        const volScalarField& scField = source().getScalarField(scName);

        word mscName = mirrorName(scName);

        volScalarField& mscField = source().getScalarField(mscName);

        volScalarField dsc(scField - mscField);

        if (uList.size() > 0)
        {
            for (unsigned i = 0; i < uList.size(); ++i)
            {
                label ic = uList[i];

                scalar msc;
                if (!Pstream::parRun())
                {
                    msc = scField.primitiveField()[ic];
                }
                else
                {
                    msc = sourceScalarFields_[scName][ic];
                }

                mscField.primitiveFieldRef()[i] = msc;
            }
        }

        if (uBList.size() > 0)
        {
            Info<< "Search boundary field..." << endl;

            std::map<word, std::vector<label>>::iterator it;
            for (it = uBList.begin(); it != uBList.end(); it++)
            {
                word pName = it->first;
                word pKey = pName + ',' + scName;    // parallel

                label ir = patchId(pName);
                if (ir < 0)
                {
                    Info<< "Error occurred in boundary mapping." << endl;

                    continue;
                }

                const std::vector<label>& bList = it->second;

                for (unsigned j = 0; j < bList.size(); j++)
                {
                    label jc = bList[j];

                    scalar sc;
                    if (!Pstream::parRun())
                    {
                        sc = scField.boundaryField()[ir][jc];
                    }
                    else
                    {
                        sc = bSourceScalarFields_[pKey][jc];
                    }

                    mscField.boundaryFieldRef()[ir][j] = sc;
                }
            }
        }
    }
}


// ************************************************************************* //
