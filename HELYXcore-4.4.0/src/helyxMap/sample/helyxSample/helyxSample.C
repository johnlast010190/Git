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

#include <dirent.h>
#include "sample/helyxSample/helyxSample.H"
#include "global/argList/argList.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "fields/fvsPatchFields/fvsPatchField/fvsPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::word Foam::helyxSample::databaseDir_;


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::getFiles(wordList& v, const word& root, const word& filter)
{
    DIR* dirp = opendir(root.c_str());

    struct dirent* dp;
    struct stat s;

    while ((dp = readdir(dirp)) != nullptr)
    {
        word dname = dp->d_name;
        if (dname == "." || dname == "..") continue;

        if (stat(dname.c_str(), &s) == 0)
        {
            if (s.st_mode)
            {
                word path = root + dname;
                if (!isFile(path.c_str())) continue;

                if (filter == "none")
                {
                    v.append(dname);
                }
                else
                {
                    label id = dname.find(filter);
                    if (id >= 0)
                    {
                        v.append(dname);
                    }
                }
            }
        }
    }

    closedir(dirp);
}


void Foam::readDirectory(wordList& v, const word& root, const word& filter)
{
    DIR* dirp = opendir(root.c_str());

    struct dirent* dp;
    struct stat s;

    while ((dp = readdir(dirp)) != nullptr)
    {
        word dname = dp->d_name;
        if (dname == "." || dname == "..") continue;

        if (stat(dname.c_str(), &s) == 0)
        {
            if (s.st_mode & S_IFDIR)
            {
                if (filter == "none")
                {
                    v.append(dname);
                }
                else
                {
                    label id = dname.find(filter);
                    if (id >= 0)
                    {
                        v.append(dname);
                    }
                }
             }
        }
    }

    closedir(dirp);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::helyxSample::writeBanner
(
    Ostream& os,
    const word classType,
    const word objectName
)
{
    os  << "/*--------------------------------*- C++ -*----------------------------------*\\" << endl;
    os  << "|       o        |                                                            |" << endl;
    os  << "|    o     o     |  HELYX (R) : Open-source CFD for Enterprise                |" << endl;
    os  << "|   o   O   o    |  Version : Dev                                             |" << endl;
    os  << "|    o     o     |  ENGYS Ltd. <http://engys.com/>                            |" << endl;
    os  << "|       o        |                                                            |" << endl;
    os  << "\\*---------------------------------------------------------------------------*/" << endl;
    os  << IOobject::foamFile << "\n{\n";
    os  << "version     2.0;\n";
    os  << "format      ascii;\n";
    os  << "class       " << classType << ";\n";
    os  << "object      " << objectName << ";\n}\n";
    os  << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::helyxSample::init()
{
    nx_ = 100;
    ny_ = 100;
    nz_ = 1;

    toDatabase_ = false;
    visualize_ = true;
    writeSolidMask_ = false;
    writeSamples_ = true;

    bodySurfaceType_ = "bodySurface";

    mapBoundary_ = true;

    mapCase_ = cwd();
    sourceCase_ = cwd();

    alphaMax_ = 0.5;

    byPhysicalType_ = false;
    byName_ = false;
    byType_ = false;
    byNameList_ = true;

    nProc_ = 1;

    word homedir = getenv("HOME");

    databaseDir_ = homedir + "/helyxAI/database";
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::helyxSample::helyxSample()
:
    helyxMap()
{
    init();
}


Foam::helyxSample::helyxSample(const Time* runTime, const fvMesh* mesh)
:
    helyxMap(),
    mesh_(mesh),
    runTime_(runTime)
{
    init();
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::helyxSample::~helyxSample()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::helyxSample::buildSearchTrees()
{
    source().constructInternalKnn(-0.1);

    farray2d trainData;

    forAll(mesh_->boundaryMesh(), r)
    {
        word type = mesh_->boundaryMesh()[r].type();
        if (type == "empty" || type == "symmetryPlane") continue;

        word pType = mesh_->boundaryMesh()[r].physicalType();
        word pName = mesh_->boundaryMesh()[r].name();

        if (!byPhysicalType_)
        {
            pType = mesh_->boundaryMesh()[r].type();
        }

        // Physical type as an indicator of body surface
        if (byType_ || byPhysicalType_)
        {
            if (pType != bodySurfaceType_) continue;
        }
        else if (byName_)
        {
            word key = "none";
            if (bodyPatchNames_.size() >= 1)
            {
                key = bodyPatchNames_[0];
            }

            if (pName != key) continue;
        }
        else if (byNameList_)
        {
            bool found = false;

            for (label ip = 0; ip < bodyPatchNames_.size(); ip++)
            {
                word filter = bodyPatchNames_[ip];
                if (pName == filter)
                {
                    found = true;
                    break;
                }

                label id = pName.find(filter);
                if (id >= 0)
                {
                    found = true;
                    break;
                }
            }

            if (!found) continue;
        }

        forAll(mesh_->boundaryMesh()[r], j)
        {
            point pt = mesh_->Cf().boundaryField()[r][j];

            scalar xb = source().xbar(pt.x());
            scalar yb = source().ybar(pt.y());
            scalar zb = source().zbar(pt.z());

            std::vector<scalar> pts(3);
            pts[0] = xb;
            pts[1] = yb;
            pts[2] = zb;

            std::pair<label, label> faceid(r, j);
            bodyFaces_.append(faceid);

            trainData.push_back(pts);
        }
    }

    if (trainData.size() < 2 && !parRun())
    {
        FatalErrorInFunction
            << "Error constructing solid object search tree, "
            << "zero number of solid object boundary faces obtained." << nl
            << "Please check whether 'bodySurfaceType' name is set correctly."
            << exit(FatalError);
    }

    if (trainData.size() >= 1)
    {
        Info<< "Body surface faces: " << trainData.size() << endl;

        bodyBoundaryTree_.build(trainData);
    }
    else
    {
        bodyBoundaryTree_.deactivate();
    }
}


void Foam::helyxSample::getGeomSample(gridField& grid)
{
    grid.setFieldSize();

    for (label i = 0; i < grid.nx_ + 1; ++i)
    {
        for (label j = 0; j < grid.ny_ + 1; ++j)
        {
            for (label k = 0; k < grid.nzpt(); ++k)
            {
                point pt = grid.gridNodes_(i, j, k);

                scalar xb = source().xbar(pt.x());
                scalar yb = source().ybar(pt.y());
                scalar zb = source().zbar(pt.z());

                std::vector<scalar> qvect(3);
                qvect[0] = xb;
                qvect[1] = yb;
                qvect[2] = zb;

                std::vector<label> nb(2);
                std::vector<scalar> dist(2);

                source().internalMap_.search(qvect, nb, dist);

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
                    inSolid = 0;
                }
                else
                {
                    if (dist1[0] > dist[0])
                    {
                        inSolid = 0;    // Nearer to fluid cell
                    }
                    else
                    {
                        // dist1[1] <= dist[0] can also be a fluid cell
                        // if cell_alpha > 1.0e-3
                        label bfaceid = nb1[0];

                        label r = bodyFaces_[bfaceid].first;
                        label jf = bodyFaces_[bfaceid].second;

                        point pf = mesh_->Cf().boundaryField()[r][jf];

                        label ic = nb[0];
                        point pc = mesh_->C()[ic];
                        scalar delta = (pc - pt) & (pf - pt);

                        if (delta > 0)
                        {
                            inSolid = 1;
                        }
                        else
                        {
                            inSolid = 0;    // delta < 0, in fluid
                        }
                    }
                }

                grid.inSolids_(i, j, k) = inSolid;
            }
        }
    }
}


void Foam::helyxSample::searchNodeField(gridField& grid, const word& funct)
{
    // Grid map does not need wall distance

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

                scalar xb = source().xbar(pt.x());
                scalar yb = source().ybar(pt.y());
                scalar zb = source().zbar(pt.z());

                std::vector<scalar> qvect(3);
                qvect[0] = xb;
                qvect[1] = yb;
                qvect[2] = zb;

                std::vector<label> nb(2);
                std::vector<scalar> dist(2);

                source().internalMap_.search(qvect, nb, dist);

                // nb stores the nearest cell index.
                // Search for the nearest body surface face center.
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
                    inSolid = 0;    // cannot find body boundary match
                }
                else
                {
                    if (dist1[0] > dist[0])
                    {
                        inSolid = 0;    // nearer to fluid cell
                        cnt1++;
                    }
                    else
                    {
                        // dist1[1] <= dist[0] can also be a fluid cell
                        // if cell_alpha > 1.0e-3
                        label bfaceid = nb1[0];

                        label r = bodyFaces_[bfaceid].first;
                        label jf = bodyFaces_[bfaceid].second;

                        point pf = mesh_->Cf().boundaryField()[r][jf];

                        label ic = nb[0];
                        point pc = mesh_->C()[ic];
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

    Info<< "Number of solid node: "
        << cnt << " " << cnt1 << " " << cnt2 << endl;
}


void Foam::helyxSample::getSampleData(const word& funct)
{
    std::map<word, gridField>::iterator it;

    for (it = sampleGrids_.begin(); it != sampleGrids_.end(); it++)
    {
        searchNodeField(it->second, funct);
    }
}


void Foam::helyxSample::setOptions(const argList& args)
{
    helyxMap::setOptions(args);

    args.optionReadIfPresent("nx", nx_);
    args.optionReadIfPresent("ny", ny_);
    args.optionReadIfPresent("nz", nz_);

    args.optionReadIfPresent("bodySurfaceType", bodySurfaceType_);
    args.optionReadIfPresent("bodySurfaceName", bodySurfaceType_);
}


void Foam::helyxSample::getBodyPatchNames()
{
    forAll(mesh_->boundaryMesh(), r)
    {
        word pType = mesh_->boundaryMesh()[r].physicalType();
        word pName = mesh_->boundaryMesh()[r].name();

        if (!byPhysicalType_)
        {
            pType = mesh_->boundaryMesh()[r].type();
        }

        if (pType == bodySurfaceType_)
        {
            bodyPatchNames_.append(pName);
        }
    }
}


void Foam::helyxSample::writeVectorField(gridField& grid, label k)
{
    if (!visualize_)
    {
        return;
    }

    word casedir = cwd();

    word visudir = casedir + "/visualization/";
    word subdir = visudir + grid.name() + "_" + runTime_->timeName();
    word runtimdir = subdir + "/" + runTime_->timeName();

    for (label iv = 0; iv < grid.vectorSamples_.size(); iv++)
    {
        word name = grid.vectorSamples_[iv];
        word fname = runtimdir + "/" + name;

        OFstream os(fname.c_str());
        if (!os)
        {
            FatalErrorInFunction
                << "Error writing vector sample field, "
                << "file: \n" << fname << " cannot be opened."
                << exit(FatalError);
        }

        scalar coef = source().scale(name);

        writeBanner(os, "volVectorField", name);

        os  << "dimensions      [0 1 -1 0 0 0 0];\n\n\n";
        os  << "internalField   nonuniform List<vector>\n";
        os  << (grid.nx_ + 1)*(grid.ny_ + 1)*grid.nzpt() << endl;
        os  << "(\n";

        DynamicList<point> uvw;
        uvw.setSize((grid.nx_ + 1)*(grid.ny_ + 1)*grid.nzpt());

        for (label i = 0; i < grid.nx_ + 1; i++)
        {
            for (label j = 0; j < grid.ny_ + 1; j++)
            {
                for (label s = 0; s < grid.nzpt(); s++)
                {
                    point v = grid.gridVectors_(i, j, s)[iv]*coef;
                    label index = grid.ijk(i, j, s);

                    uvw[index] = v;
                }
            }
        }

        for (label i = 0; i < uvw.size(); i++)
        {
            point v = uvw[i];
            os  << "(" << v.x() << " " << v.y() << " " << v.z() << ")\n";
        }

        os  << ")\n;\n";
        os  << "boundaryField\n{\n";
        os  << "    left\n    {\n";
        os  << "    type    inlet;\n";
        os  << "    value           nonuniform List<vector>\n";
        os  << (grid.ny_ + 1)*grid.nzpt() << endl;
        os  << "(\n";

        for (label j = 0; j < grid.ny_ + 1; j++)
        {
            for (label s = 0; s < grid.nzpt(); s++)
            {
                label idx = grid.ijk(0, j, s);
                point v = uvw[idx];

                os  << "(" << v.x() << " " << v.y() << " " << v.z() << ")\n";
            }
        }

        os  << ")\n;\n}\n";
        os  << "    right\n    {\n";
        os  << "    type    outlet;\n";
        os  << "    value  nonuniform List<vector>\n";
        os  << (grid.ny_ + 1)*grid.nzpt() << endl;
        os  << "(\n";

        for (label j = 0; j < grid.ny_ + 1; j++)
        {
            for (label s = 0; s < grid.nzpt(); s++)
            {
                label idx = grid.ijk(grid.nx_, j, s);
                point v = uvw[idx];

                os  << "(" << v.x() << " " << v.y() << " " << v.z() << ")\n";
            }
        }

        os  << ")\n;\n}\n";
        os  << "    top\n    {\n";
        os  << "    type    wall;\n";
        os  << "    value  nonuniform List<vector>\n";
        os  << (grid.nx_ + 1)*grid.nzpt() << endl;
        os  <<"(\n";

        for (label i = 0; i < grid.nx_ + 1; i++)
        {
            for (label s = 0; s < grid.nzpt(); s++)
            {
                label idx = grid.ijk(i, 0, s);
                point v = uvw[idx];

                os  << "(" << v.x() << " " << v.y() << " " << v.z() << ")\n";
            }
        }

        os  << ")\n;\n}\n";
        os  << "    bottom\n    {\n";
        os  << "    type    wall;\n";
        os  << "    value   nonuniform List<vector>\n";
        os  << (grid.nx_ + 1)*grid.nzpt() << endl;
        os  << "(\n";

        for (label i = 0; i < grid.nx_ + 1; i++)
        {
            for (label s = 0; s < grid.nzpt(); s++)
            {
                label idx = grid.ijk(i, grid.ny_, s);
                point v = uvw[idx];

                os  << "(" << v.x() << " " << v.y() << " " << v.z() << ")\n";
            }
        }

        os  << ")\n;\n}\n";
        os  << "    front\n{\n";
        os  << "    type    empty;\n}\n\n";
        os  << "    back\n{\n";
        os  << "    type    empty;\n}\n}\n";
    }
}


void Foam::helyxSample::writeScalarField(gridField& grid, label k)
{
    if (!visualize_)
    {
        return;
    }

    word casedir = cwd();

    word visudir = casedir + "/visualization/";
    word subdir = visudir + grid.name() + "_" + runTime_->timeName();
    word runtimdir = subdir + "/" + runTime_->timeName();

    for (label isc = 0; isc < grid.scalarSamples_.size(); isc++)
    {
        word name = grid.scalarSamples_[isc];
        word fname = runtimdir + "/" + name;

        OFstream os(fname.c_str());
        if (!os)
        {
            FatalErrorInFunction
                << "Error writing scalar sample field, "
                << "file: \n" << fname << " cannot be opened."
                << exit(FatalError);
        }

        scalar coef = source().scale(name);

        writeBanner(os, "volScalarField", name);

        os  << "dimensions      [0 0 0 0 0 0 0];\n\n\n";
        os  << "internalField   nonuniform List<scalar>\n";
        os  << (grid.nx_ + 1)*(grid.ny_ + 1)*grid.nzpt() << endl;
        os  << "(\n";

        DynamicList<scalar> scfield;
        scfield.setSize((grid.nx_ + 1)*(grid.ny_ + 1)*grid.nzpt());

        for (label i = 0; i < grid.nx_ + 1; i++)
        {
            for (label j = 0; j < grid.ny_ + 1; j++)
            {
                for (label s = 0; s < grid.nzpt(); s++)
                {
                    scalar d = grid.gridFields_(i, j, s)[isc]*coef;
                    label index = grid.ijk(i, j, s);

                    scfield[index] = d;
                }
            }
        }

        for (label i = 0; i < scfield.size(); i++)
        {
            os  << scfield[i] << endl;
        }

        os  <<")\n;\n";
        os  <<"boundaryField\n{\n";
        os  <<"    left\n    {\n";
        os  <<"        type     zeroGradient;\n        }\n";
        os  <<"    right\n    {\n";
        os  <<"        type    zeroGradient;\n        }\n\n";

        os  <<"    top\n    {\n";
        os  <<"        type    zeroGradient;\n        }\n\n";

        os  <<"    bottom\n    {\n";
        os  <<"        type    zeroGradient;\n           }\n\n";

        os  <<"    front\n{\n";
        os  <<"    type    empty;\n}\n\n";

        os  <<"    back\n{\n";
        os  <<"    type    empty;\n}\n}\n";
    }
}


void Foam::helyxSample::writeSolidMask(gridField& grid)
{
    if (!visualize_)
    {
        return;
    }

    word casedir = cwd();

    word visudir = casedir + "/visualization/";
    word subdir = visudir + grid.name() + "_" + runTime_->timeName();
    word runtimdir = subdir + "/" + runTime_->timeName();

    word name = "solidMask";
    word fname = runtimdir + "/" + name;

    OFstream os(fname.c_str());
    if (!os)
    {
        FatalErrorInFunction
            << "Error writing scalar sample field, "
            << "file: \n" << fname << " cannot be opened."
            << exit(FatalError);
    }

    writeBanner(os, "volScalarField", name);

    os  << "dimensions      [0 0 0 0 0 0 0];\n\n\n";
    os  << "internalField   nonuniform List<scalar>\n";
    os  << (grid.nx_ + 1)*(grid.ny_ + 1)*grid.nzpt() << endl;
    os  << "(\n";

    DynamicList<scalar> scfield;
    scfield.setSize((grid.nx_ + 1)*(grid.ny_ + 1)*grid.nzpt());

    for (label i = 0; i < grid.nx_ + 1; i++)
    {
        for (label j = 0; j < grid.ny_ + 1; j++)
        {
            for (label s = 0; s < grid.nzpt(); s++)
            {
                scalar d = grid.inSolids_(i, j, s);
                label index = grid.ijk(i, j, s);

                scfield[index] = d;
            }
        }
    }

    for (label i = 0; i < scfield.size(); i++)
    {
        os  << scfield[i] << endl;
    }

    os  <<")\n;\n";
    os  <<"boundaryField\n{\n";
    os  <<"    left\n    {\n";
    os  <<"        type     zeroGradient;\n        }\n";
    os  <<"    right\n    {\n";
    os  <<"        type    zeroGradient;\n        }\n\n";

    os  <<"    top\n    {\n";
    os  <<"        type    zeroGradient;\n        }\n\n";

    os  <<"    bottom\n    {\n";
    os  <<"        type    zeroGradient;\n           }\n\n";

    os  <<"    front\n{\n";
    os  <<"    type    empty;\n}\n\n";

    os  <<"    back\n{\n";
    os  <<"    type    empty;\n}\n}\n";
}


void Foam::helyxSample::writeMeshDict(gridField& grid)
{
    if (!visualize_)
    {
        return;
    }

    scalar coef = 1.0;

    scalar xmin = grid.xmin() - 0.5*grid.dx_*coef;
    scalar xmax = grid.xmax() + 0.5*grid.dx_*coef;

    scalar ymin = grid.ymin() - 0.5*grid.dy_*coef;
    scalar ymax = grid.ymax() + 0.5*grid.dy_*coef;

    scalar zmin = grid.zmin() - 0.5*grid.dz_*coef;
    scalar zmax = grid.zmax() + 0.5*grid.dz_*coef;

    word casedir = cwd();

    word visudir = casedir + "/visualization/";
    word subdir = visudir + grid.name() + "_" + runTime_->timeName();
    word sysdir = subdir + "/system";
    word gfile = sysdir + "/blockMeshDict";

    OFstream os(gfile.c_str());
    if (!os)
    {
        FatalErrorInFunction
            << "Error writing 'blockMeshDict' file: \n"
            << gfile << " cannot be opened."
            << exit(FatalError);
    }

    writeBanner(os, "dictionary", "blockMeshDict");

    DynamicList<point> points;
    points.resize(8);

    points[0] = point(xmin, ymin, zmin);
    points[1] = point(xmax, ymin, zmin);
    points[2] = point(xmax, ymax, zmin);
    points[3] = point(xmin, ymax, zmin);
    points[4] = point(xmin, ymin, zmax);
    points[5] = point(xmax, ymin, zmax);
    points[6] = point(xmax, ymax, zmax);
    points[7] = point(xmin, ymax, zmax);

    os  << "convertToMeters 1;\n\n";
    os  << " vertices\n\(\n";

    for (label i = 0; i < points.size(); i++)
    {
        const point& pt = points[i];

        os  << "      (" << pt.x() << " " << pt.y() << " " << pt.z() << ")\n";
    }

    os  << ");\n\n";
    os  << "blocks\n\(\n";
    os  << "  hex (0 1 2 3 4 5 6 7) (" << grid.nx_ + 1 << " " << grid.ny_ + 1
        << " " << grid.nzpt() << ") simpleGrading (1 1 1)\n";
    os  << ");\n" << endl;
    os  << "edges\n(\n);\n";

    os  << "boundary\n\(\n\
        top\n\
        {\n\
            type patch;\n\
            faces\n\
            (\n\
            (3 7 6 2)\n\
            );\n\
        }\n\n";

    os  << "left\n\
        {\n\
            type patch;\n\
            faces\n\
            (\n\
                (0 4 7 3  )\n\
            );\n\
        }\n\n";

    os  << "bottom\n\
        {\n\
            type patch;\n\
            faces\n\
            (\n\
                (1 5 4 0)\n\
            );\n\
        }\n\n";

    os  << "right\n\
        {\n\
            type patch;\n\
            faces\n\
            (\n\
                (2 6 5 1)\n\
            );\n\
        }\n\n";

    os  << "front\n\
        {\n\
            type empty;\n\
            faces\n\
            (\n\
                (4 5 6 7)\n\
            );\n\
        }\n\n";

    os  << "back\n\
        {\n\
        type empty;\n\
        faces\n\
        (\n\
          (0 3 2 1)\n\
        );\n\
        }\n\
        \n\
        );\n";
}


void Foam::helyxSample::writeVisualization()
{
    if (!visualize_)
    {
        Info<< "Visualization option turned off, "
            << "visualization files will not be written." << endl;

        return;
    }

    Info<< "Writing visualization file..." << endl;

    const fileName path("visualization");
    mkDir(path);

    std::map<word, gridField>::iterator it;

    word caseDir = cwd();
    rootCase_ = caseDir;
    word visuDir = caseDir + "/visualization/";

    for (it = sampleGrids_.begin(); it != sampleGrids_.end(); it++)
    {
        gridField& grid = it->second;
        const fileName subDir
        (
            visuDir + grid.name() + "_" + runTime_->timeName()
        );
        mkDir(subDir);

        const fileName sysDir(subDir + "/system");
        mkDir(sysDir);

        const fileName constDir(subDir + "/constant");
        mkDir(constDir);

        const fileName meshDir(constDir + "/polyMesh");
        mkDir(meshDir);

        const fileName runtimeDir(subDir + "/" + runTime_->timeName());
        mkDir(runtimeDir);

        writeMeshDict(it->second);

        cp("system/controlDict", sysDir);
        cp("system/fvSolution", sysDir);
        cp("system/fvSchemes", sysDir);

        writeScalarField(it->second);

        if (writeSolidMask_)
        {
            writeSolidMask(it->second);
        }

        writeVectorField(it->second);

        chDir(subDir);

        const word cmd = "blockMesh";
        system(cmd);

        chDir(rootCase_);
    }
}


void Foam::helyxSample::writeNodeField()
{
    if (!writeSamples_)
    {
        Info<< "Write samples option set to false. "
            << "Sample data will not be written into file." << endl;

        return;
    }

    std::map<word, gridField>::iterator it;
    for (it = sampleGrids_.begin(); it != sampleGrids_.end(); it++)
    {
        writeNodeField(it->second);
    }
}


void Foam::helyxSample::saveToDatabase(const word& funct)
{
    if (!toDatabase_)
    {
        Info<< "Option to_database is false, "
            << "data will not be saved to database." << endl;

        return;
    }

    std::map<word, gridField>::iterator it;
    for (it = sampleGrids_.begin(); it != sampleGrids_.end(); it++)
    {
        saveToDatabase(it->second, funct);
    }
}


void Foam::helyxSample::saveObjectMasks(const word& type)
{
    std::map<word, gridField>::iterator it;

    for (it = sampleGrids_.begin(); it != sampleGrids_.end(); it++)
    {
        saveObjectMask(it->second, type);
    }
}


void Foam::helyxSample::saveObjectMask(gridField& grid, const word& type)
{
    word casedir = cwd();

    word dbsdir = databaseDir_ + "/" + type;
    word dbsfile = dbsdir + "/objectMasks.dbs";
    databasePath_ = dbsfile;

    word fcaselist = dbsdir + "/caselist.dbs";
    std::vector<word> cases;

    std::ifstream is(fcaselist.c_str());
    if (is)
    {
        word str;

        while (!is.eof())
        {
            std::getline(is, str, '\n');
            if (str.length() < 2)
            {
                continue;
            }

            cases.push_back(str);
        }

        is.close();
    }

    bool exists = false;
    if (cases.size() >= 1)
    {
        std::vector<word>::iterator it;
        it = std::find(cases.begin(), cases.end(), casedir);

        if (it != cases.end())
        {
            exists = true;
        }
    }

    if (exists)
    {
        Info<< "Warning: current case already in database." << endl;

        return;
    }

    OFstream os1
    (
        fcaselist.c_str(),
        IOstream::ASCII,
        IOstream::currentVersion,
        runTime_->writeCompression(),
        true
    );

    os1 << casedir << endl;

    const fileName path(dbsdir);
    mkDir(path);

    OFstream os
    (
        dbsfile.c_str(),
        IOstream::ASCII,
        IOstream::currentVersion,
        runTime_->writeCompression(),
        true
    );

    word desc =
        type + "," + casedir + "," + grid.name() +
        "," + std::to_string(grid.nx_+1) +
        "," + std::to_string(grid.ny_+1) +
        "," + std::to_string(grid.nzpt()) +
        "," + std::to_string(grid.xmin()) +
        "," + std::to_string(grid.ymin()) +
        "," + std::to_string(grid.zmin()) +
        "," + std::to_string(grid.xmax()) +
        "," + std::to_string(grid.ymax()) +
        "," + std::to_string(grid.zmax());

    os  << desc << endl;

    word cnts = "mask";

    for (label i = 0; i < grid.nx_ + 1; ++i)
    {
        for (label j = 0; j < grid.ny_ + 1; ++j)
        {
            for (label k = 0; k < grid.nzpt(); ++k)
            {
                label inSolid = grid.inSolids_(i, j, k);

                cnts += ("," + std::to_string(inSolid));
            }
        }
    }

    os  << cnts << endl;
}


void Foam::helyxSample::saveToDatabase(gridField& grid, const word& funct)
{
    if (parRun())
    {
        return;
    }

    word caseFileList = databaseDir_ + "/caselist.dbs";

    // Check whether case name is already in database
    word bkName = caseName() + ":" + grid.name();

    if (funct == "geometry")
    {
        // Save case geometry
        Info<< "Write geometric information into database." << endl;
        word geofile = databaseDir_ + "/geom.dbs";

        OFstream os
        (
            geofile.c_str(),
            IOstream::ASCII,
            IOstream::currentVersion,
            runTime_->writeCompression(),
            true
        );

        if (!os)
        {
            FatalErrorInFunction
                << "Error saving sample geometry to database, "
                << "database file \n" << geofile << " cannot be opened."
                << exit(FatalError);
        }

        os  << bkName << endl;
        os  << grid.size() << " " << grid.nx_ << " "
            << grid.ny_ << " " << grid.nz_ << endl;

        std::vector<label> inSolid(grid.size());

        label cnt = 0;
        for (label i = 0; i < grid.nx_ + 1; ++i)
        {
            for (label j = 0; j < grid.ny_ + 1; ++j)
            {
                for (label k = 0; k < grid.nzpt(); ++k)
                {
                    label idx = grid.ijk(i, j, k);

                    inSolid[idx] = grid.inSolids_(i, j, k);

                    cnt += inSolid[idx];
                }
            }
        }

        Info<< "In solids: " << cnt << endl;

        for (unsigned i = 0; i < inSolid.size(); i++)
        {
            os  << inSolid[i] << " ";
        }
        os  << endl;

        return;
    }

    std::ifstream is(caseFileList.c_str());
    std::map<word, word> savedCases;

    if (is)
    {
        word str;

        while (!is.eof())
        {
            std::getline(is, str, '\n');
            if (str.length() < 3) continue;

            wordList vtmp;
            stringTokenize(vtmp, str, "=");

            if (vtmp.size() != 2) continue;

            savedCases[vtmp[0]] = vtmp[1];
        }

        is.close();
    }

    bool caseExists = true;
    if (savedCases.size() == 0)
    {
        caseExists = false;
    }
    else
    {
        std::map<word, word>::iterator it = savedCases.find(bkName);
        if (it == savedCases.end())
        {
            caseExists = false;
        }
    }

    if (!caseExists)
    {
        OFstream os
        (
            caseFileList.c_str(),
            IOstream::ASCII,
            IOstream::currentVersion,
            runTime_->writeCompression(),
            true
        );

        if (!os)
        {
            FatalErrorInFunction
                << "Error saving sample data to database, "
                << "database file \n" << caseFileList << " cannot be opened."
                << exit(FatalError);
        }

        os  << bkName << "=" << source().rhoRef_ << " " << source().Uref_.x()
            << " " << source().Uref_.y() << " " << source().Uref_.z() << " "
            << grid.xmin() << " " << grid.ymin() << " " <<grid.zmin() << " "
            << grid.xmax() << " " << grid.ymax() << " " <<grid.zmax() << endl;
    }

    // Write field data
    const fileName caseDatadir(databaseDir_ + "/caseData/");
    mkDir(caseDatadir);

    word datafile = caseDatadir + bkName + ".dbs";

    OFstream os(datafile.c_str());
    if (!os)
    {
        FatalErrorInFunction
            << "Error saving sample data to database, "
            << "database file \n" << datafile << " cannot be opened."
            << exit(FatalError);
    }

    os  << grid.nx_ << " " << grid.ny_ << " " << grid.nz_ << " "
        << grid.containSolidObject_ << " " << grid.xmin() << " "
        << grid.ymin() << " " << grid.zmin() << " " << grid.xmax() << " "
        << grid.ymax() << " " << grid.zmax() << endl;

    word header = "insolid ";

    forAllConstIter(HashSet<word>, mapScalarFields_, iter)
    {
        header += (iter() + " ");
    }
    forAllConstIter(HashSet<word>, mapVectorFields_, iter)
    {
        const word& vname = iter();

        word vx = vname + "_x";
        word vy = vname + "_y";
        word vz = vname + "_z";

        header += (" " + vx + " " + vy + " " + vz);
    }

    os  << header << endl;

    DynamicList<std::vector<scalar>> data;
    label sz = grid.size();
    data.resize(sz);

    label ndata = 1 + mapScalarFields_.size() + 3*mapVectorFields_.size();

    for (label i = 0; i < grid.nx_ + 1; i++)
    {
        for (label j = 0; j < grid.ny_ + 1; j++)
        {
            for (label k = 0; k < grid.nzpt(); k++)
            {
                label idx = grid.ijk(i, j, k);

                std::vector<scalar> vect(ndata);
                vect[0] = grid.inSolids_(i, j, k);

                for (label isc = 0; isc < mapScalarFields_.size(); isc++)
                {
                    scalar d = grid.gridFields_(i, j, k)[isc];
                    vect[isc + 1] = d;
                }

                if (mapVectorFields_.size() >= 1)
                {
                    for (label iv = 0; iv < mapVectorFields_.size(); iv++)
                    {
                        point vd = grid.gridVectors_(i, j, k)[iv];
                        label idata = mapScalarFields_.size() + 3*iv + 1;

                        vect[idata] = vd.x();
                        vect[idata + 1] = vd.y();
                        vect[idata + 2] = vd.z();
                    }
                }

                data[idx] = vect;
            }
        }
    }

    os  << data.size() << endl;

    os.precision(4);

    for (label ik = 0; ik < data.size(); ik++)
    {
        for (label id = 0; id < ndata; id++)
        {
            os  << data[ik][id];

            if (id < ndata - 1)
            {
                os  << " ";
            }
            else
            {
                os  << endl;
            }
        }
    }
}


void Foam::helyxSample::writeNodeField(gridField& grid)
{
    Info<< "mapScalarFields: " << mapScalarFields_.size() << endl;

    const fileName dir(cwd() + "/postProcessing/");

    word filename = dir + grid.name() + "_" + runTime_->timeName();
    Info<< filename << endl;

    mkDir(dir);

    label sz = grid.size();

    // Write the distance file
    word dfile = filename + "_dists";
    OFstream osd(dfile.c_str());

    std::vector<word> ddata(sz);

    if (osd)
    {
        word ss;
        for (label i = 0; i < grid.nx_ + 1; i++)
        {
            for (label j = 0; j < grid.ny_ + 1; j++)
            {
                for (label k = 0; k < grid.nzpt(); k++)
                {
                    label idx = grid.ijk(i, j, k);
                    ddata[idx] = doubleToString(grid.dists_(i, j, k));
                }
            }
        }

        for (label s = 0; s < sz; s++)
        {
            ss += (ddata[s] + "\n");
        }

        osd << ss;
    }

    OFstream os(filename.c_str());
    if (!os)
    {
        FatalErrorInFunction
            << "Error saving sample data, "
            << "sample file \n" << filename << " cannot be opened."
            << exit(FatalError);
    }

    word header = "x y z";

    forAllConstIter(HashSet<word>, mapScalarFields_, iter)
    {
        header += (" " + iter());
    }
    forAllConstIter(HashSet<word>, mapVectorFields_, iter)
    {
        const word& vname = iter();

        word vx = vname + "_x";
        word vy = vname + "_y";
        word vz = vname + "_z";

        header += (" " + vx + " " + vy + " " + vz);
    }

    header += " in_solid_object";
    Info<< header << endl;

    os  << header << endl;

    DynamicList<std::vector<scalar>> data;
    data.resize(sz);

    label ndata = 3 + mapScalarFields_.size() + 1 + 3*mapVectorFields_.size();

    for (label i = 0; i < grid.nx_ + 1; i++)
    {
        for (label j = 0; j < grid.ny_ + 1; j++)
        {
            for (label k = 0; k < grid.nzpt(); k++)
            {
                label idx = grid.ijk(i, j, k);
                std::vector<scalar> vect(ndata);

                vect[0] = grid.gridNodes_(i, j, k).x();
                vect[1] = grid.gridNodes_(i, j, k).y();
                vect[2] = grid.gridNodes_(i, j, k).z();

                forAllConstIter(HashSet<word>, mapScalarFields_, iter)
                {
                    const word& name = iter();

                    label isc = grid.scalarId(name);
                    scalar coef = source().scale(name);
                    scalar d = grid.gridFields_(i, j, k)[isc];

                    vect[3 + isc] = d*coef;
                }

                forAllConstIter(HashSet<word>, mapVectorFields_, iter)
                {
                    const word& name = iter();

                    label iv = grid.vectorId(name);
                    point vd = grid.gridVectors_(i, j, k)[iv];
                    scalar coef = source().scale(name);

                    vd *= coef;

                    label idata = 3 + mapScalarFields_.size() + 3*iv;

                    vect[idata] = vd.x();
                    vect[idata + 1] = vd.y();
                    vect[idata + 2] = vd.z();
                }

                vect[ndata - 1] = grid.inSolids_(i, j, k);

                data[idx] = vect;
            }
        }
    }

    os.precision(4);

    for (label ik = 0; ik < data.size(); ik++)
    {
        for (label id = 0; id < ndata; id++)
        {
            os  << data[ik][id];

            if (id < ndata - 1)
            {
                os  << " ";
            }
            else
            {
                os  << endl;
            }
        }
    }
}


void Foam::helyxSample::setInput(const dictionary& dict, const argList& args)
{
    helyxMap::setInput(dict, args);

    getFieldTypes(*mesh_, runTime_->timeName());

    visualize_ = dict.lookupOrDefault<bool>("visualize", false);

    writeSolidMask_ =
        visualize_
      ? dict.lookupOrDefault<bool>("writeSolidMask", false)
      : false;

    bodySurfaceType_ = dict.lookupOrDefault<word>("bodySurfaceType", "car");

    toDatabase_ = dict.lookupOrDefault<bool>("to_database", false);

    word identifyObject =
        dict.lookupOrDefault<word>("identifyObject", "byNameList");

    if (identifyObject == "byType")
    {
        byPhysicalType_ = false;
        byType_ = true;
        byName_ = false;
        byNameList_ = false;
    }
    else if (identifyObject == "byPhysicalType")
    {
        byPhysicalType_ = true;
        byType_ = false;
        byName_ = false;
        byNameList_ = false;
    }
    else if (identifyObject == "byName")
    {
        byPhysicalType_ = false;
        byType_ = false;
        byName_ = true;
        byNameList_ = false;
    }
    else if (identifyObject == "byNameList")
    {
        byPhysicalType_ = false;
        byType_ = false;
        byName_ = false;
        byNameList_ = true;
    }

    word blockName = dict.lookupOrDefault<word>("blockName", "block1");

    point p0(0, 0, 0);
    point p1(1, 1, 1);

    boundBox bbox0(p0, p1);
    boundBox bbox = dict.lookupOrDefault<boundBox>("boundBox", bbox0);

    point dp = bbox.max() - bbox.min();

    scalar coef = 0.05;

    point pmin = bbox.min() - coef*dp;
    point pmax = bbox.max() + coef*dp;

    sampleBox_ = new boundBox(pmin, pmax);

    vector ndivs0(127, 127, 1);
    vector ndivisions = dict.lookupOrDefault<vector>("ndivisions", ndivs0);

    label nx = label(ndivisions[0]);
    label ny = label(ndivisions[1]);
    label nz = label(ndivisions[2]);

    bool containSolidObject = true;

    writeSamples_= dict.lookupOrDefault<bool>("write_samples", true);
    bool interp = dict.lookupOrDefault<bool>("interpolation", false);

    gridField grid(bbox, nx, ny, nz, containSolidObject);

    grid.interp_ = interp;
    grid.name_ = blockName;

    grid.setSampleNames(mapScalarFields_, mapVectorFields_);
    grid.generateNodes();
    grid.setFieldSize();

    sampleGrids_[blockName] = grid;

    const polyBoundaryMesh& bMesh = mesh_->boundaryMesh();

    labelHashSet includePatches;

    if (dict.found("patches"))
    {
        includePatches = bMesh.patchSet(wordReList(dict.lookup("patches")));
    }

    const wordList allPatchNames(bMesh.names());

    forAllConstIter(labelHashSet, includePatches, it)
    {
        bodyPatchNames_.append(allPatchNames[*it]);
    }
}


// ************************************************************************* //
