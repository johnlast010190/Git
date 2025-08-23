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

#include "sample/helyxSamplePar/helyxSamplePar.H"
#include "global/argList/argList.H"
#include "fields/fvsPatchFields/fvsPatchField/fvsPatchFields.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::interProcTrans
(
    DynamicList<T>& field,
    const DynamicList<T>& source,
    const word& name
)
{
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // All send to master
    if (Pstream::myProcNo() != Pstream::masterNo())
    {
        UOPstream toMaster(Pstream::masterNo(), pBufs);
        toMaster << source;
    }

    pBufs.finishedSends();

    if (Pstream::myProcNo() == Pstream::masterNo())
    {
        // Collect my own data
        field.append(source);

        for
        (
            int slave = Pstream::firstSlave();
            slave <= Pstream::lastSlave();
            slave++
        )
        {
            DynamicList<T> posGlob;

            UIPstream fromSlave(slave, pBufs);
            fromSlave >> posGlob;

            field.append(posGlob);
        }
    }

    Info<< "Field transferred for: " << name << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::helyxSamplePar::helyxSamplePar()
:
    helyxSample(),
    sampleMethod_("nearest")
{}


Foam::helyxSamplePar::helyxSamplePar(const Time* runTime, const fvMesh* mesh)
:
    helyxSample(runTime, mesh),
    sampleMethod_("nearest")
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::helyxSamplePar::~helyxSamplePar()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::helyxSamplePar::setOptions(const argList& args)
{
    helyxSample::setOptions(args);

    args.optionReadIfPresent("sampleMethod", sampleMethod_);
}


void Foam::helyxSamplePar::getDirNames(wordList& files, const word& pattern)
{
    getFiles(files, pattern);
}


void Foam::helyxSamplePar::getMapTime()
{
    if (mapTimeName_ != "latest")
    {
        return;
    }

    word cdir = caseDir_ + "/processor" + std::to_string(myId_);
    wordList dirs;

    word pattern = cdir + "/*";

    getDirNames(dirs, pattern);

    scalar latestTime = 0;

    for (label i = 0; i < dirs.size(); i++)
    {
        if (!isDir(dirs[i].c_str()))
        {
            continue;
        }

        wordList vtmp;
        stringTokenize(vtmp, dirs[i], "/");

        word wlast = vtmp[vtmp.size() - 1];
        label itime = atoi(wlast.c_str());

        latestTime = max(latestTime, itime);
    }

    mapTimeName_ = name(latestTime);
}


void Foam::helyxSamplePar::getBodySurface()
{
    if (!Pstream::parRun())
    {
        return;
    }

    getBodyPatchNames();

    List<vectorField> posLoc(Pstream::nProcs());
    label myid = Pstream::myProcNo();

    DynamicList<point> xyz;

    forAll(mesh_->boundaryMesh(), r)
    {
        word type = mesh_->boundaryMesh()[r].type();

        if (type == "empty" || type == "symmetryPlane")
        {
            continue;
        }

        word name = mesh_->boundaryMesh()[r].name();

        if (!inList(name, bodyPatchNames_))
        {
            continue;
        }

        forAll(mesh_->boundaryMesh()[r], j)
        {
            const point& pt = mesh_->Cf().boundaryField()[r][j];
            xyz.append(pt);
        }
    }

    posLoc[myid] = xyz;

    Pstream::allGatherList(posLoc);

    bodySurfCenter_ =
        ListListOps::combine<vectorField>(posLoc, accessOp<vectorField>());
}


void Foam::helyxSamplePar::combineFieldsBK()
{
    if (!Pstream::parRun())
    {
        return;
    }

    List<vectorField> posLoc(Pstream::nProcs());
    label myid = Pstream::myProcNo();

    DynamicList<point> xyz;

    xyz.setSize(mesh_->C().size());
    forAll(mesh_->C(), i)
    {
        xyz[i] = mesh_->C()[i];
    }

    posLoc[myid] = xyz;

    Pstream::allGatherList(posLoc);

    if (master())
    {
        sourceXyz_=
            ListListOps::combine<vectorField>(posLoc, accessOp<vectorField>());
    }

    wait();

    forAllConstIter(HashSet<word>, mapScalarFields_, iter)
    {
        const word& name = iter();

        const volScalarField& sfield = source().getScalarField(name);
        label sz = sfield.primitiveField().size();

        scalarField scfldLoc(sz);
        for (label i = 0; i < sz; i++)
        {
            scfldLoc[i] = sfield[i];
        }

        List<scalarField> sField(Pstream::nProcs());
        sField[myid] = scfldLoc;
        Pstream::allGatherList(sField);

        scalarField fldGlobal =
            ListListOps::combine<scalarField>(sField, accessOp<scalarField>());

        sourceScalarFields_[name].resize(0);
        if (master())
        {
            sourceScalarFields_[name] = fldGlobal;
        }

        wait();

        fldGlobal.resize(0);
    }

    // Velocity fields
    forAllConstIter(HashSet<word>, mapVectorFields_, iter)
    {
        const word& name = iter();

        const volVectorField& vField = source().getVectorField(name);
        label sz = vField.primitiveField().size();

        vectorField vectLoc(sz);
        for (label i = 0; i < sz; i++)
        {
            vectLoc[i] = vField[i];
        }

        List<vectorField> vecField(Pstream::nProcs());
        vecField[myid] = vectLoc;
        Pstream::allGatherList(vecField);

        vectorField vfldGlobal =
            ListListOps::combine<vectorField>
            (
                vecField,
                accessOp<vectorField>()
            );

        sourceVectorFields_[name].resize(0);
        if (master())
        {
            sourceVectorFields_[name] = vfldGlobal;
        }

        wait();

        vfldGlobal.resize(0);
    }
}


void Foam::helyxSamplePar::combineFields()
{
    if (!Pstream::parRun())
    {
        return;
    }

    // Cell centre coordinates
    interProcTrans(sourceXyz_, source().xyz(), "cellVerts");

    // Scalar fields
    forAllConstIter(HashSet<word>, mapScalarFields_, iter)
    {
        const word& name = iter();
        const DynamicList<scalar>& scfld = source().scalarFields_[name];

        interProcTrans(sourceScalarFields_[name], scfld, name);
    }

    // Vector fields
    forAllConstIter(HashSet<word>, mapVectorFields_, iter)
    {
        const word& name = iter();
        const DynamicList<point>& vecfld = source().vectorFields_[name];

        interProcTrans(sourceVectorFields_[name], vecfld, name);
    }
}


void Foam::helyxSamplePar::constructKnn()
{
    if (!master())
    {
        return;
    }

    label ncells = sourceXyz_.size();
    label cnt = 0;
    farray2d trainData;

    for (label i = 0; i < ncells; i++)
    {
        const point& pt = sourceXyz_[i];

        std::vector<scalar> pts(3);
        pts[0] = pt.x();
        pts[1] = pt.y();
        pts[2] = pt.z();

        trainData.push_back(pts);

        cellMap_[cnt] = i;

        cnt++;
    }

    Info<< "Build knn... Size: " << trainData.size() << endl;
    internalMap_.build(trainData);

    trainData.resize(0);
    for (label i = 0; i < bodySurfCenter_.size(); i++)
    {
        const point& pt = bodySurfCenter_[i];

        std::vector<scalar> pts(3);
        pts[0] = pt.x();
        pts[1] = pt.y();
        pts[2] = pt.z();

        trainData.push_back(pts);
    }

    Info<< "Build object knn... Size: " << trainData.size() << endl;
    bodyBoundaryTree_.build(trainData);
}


void Foam::helyxSamplePar::searchNodeField(gridField& grid, const word& funct)
{
    if (!master())
    {
        return;
    }

    // Set placeholders
    grid.setFieldSize();

    label cnt = 0;
    label cnt1 = 0;
    label cnt2 = 0;

    for (label i = 0; i < grid.nx_ + 1; ++i)
    {
        for (label j = 0; j < grid.ny_ + 1; ++j)
        {
            for (label k = 0; k < grid.nzpt(); ++k)
            {
                const point& pt = grid.gridNodes_(i, j, k);

                std::vector<scalar> qvect(3);
                qvect[0] = pt.x();
                qvect[1] = pt.y();
                qvect[2] = pt.z();

                std::vector<label> nb;
                std::vector<scalar> dist;

                internalMap_.search(qvect, nb, dist);

                // Search for the nearest body surface face center
                std::vector<label> nb1;
                std::vector<scalar> dist1;

                if (bodyBoundaryTree_.active)
                {
                    bodyBoundaryTree_.search(qvect, nb1, dist1);
                }

                // Check if a grid node is in solid
                label inSolid = 1;

                if (nb1.size() == 0)
                {
                    Info<< "Cannot find boundary match..." << endl;
                    inSolid = 0;
                }
                else
                {
                    if (dist1[0] > dist[0])
                    {
                        // Nearer to fluid cell
                        inSolid = 0;
                        cnt1++;
                    }
                    else
                    {
                        // dist1[1] <= dist[0] can also be a fluid cell
                        // if cell_alpha > 1.0e-3
                        label bid = nb1[0];
                        const point& pf = bodySurfCenter_[bid];

                        label ic = nb[0];
                        point pc = sourceXyz_[ic];

                        scalar delta = (pc - pt) & (pf - pt);
                        if (delta > 0)
                        {
                            inSolid = 1;
                        }
                        else
                        {
                            inSolid = 0;    // delta < 0, in fluid
                            cnt2++;
                        }
                    }
                }

                grid.inSolids_(i, j, k) = inSolid;
                cnt += inSolid;

                if (funct == "field")
                {
                    if (grid.scalarSamples_.size() >= 1)
                    {
                        getGridField(grid, i, j, k, nb, dist, inSolid);
                    }
                    if (grid.vectorSamples_.size() >= 1)
                    {
                        getGridVectorField(grid, i, j, k, nb, dist, inSolid);
                    }
                }
            }
        }
    }

    Info<< "Number of node in solids: " << cnt << " "
        << cnt1 << " " << cnt2 << endl;
}


// ************************************************************************* //
