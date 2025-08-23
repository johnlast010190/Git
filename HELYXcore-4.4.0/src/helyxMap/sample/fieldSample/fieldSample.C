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

#include "sample/fieldSample/fieldSample.H"
#include "helyxMap.H"
#include "fields/fvsPatchFields/fvsPatchField/fvsPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldSample::fieldSample()
:
    gridField(),
    bodySurfaceType_("bodySurface"),
    byName_(true),
    byPhysicalType_(false)
{
    nx_ = 100;
    ny_ = 100;
    nz_ = 1;
}


Foam::fieldSample::fieldSample
(
    label nx,
    label ny,
    label nz,
    scalar x0,
    scalar y0,
    scalar z0,
    scalar x1,
    scalar y1,
    scalar z1
)
:
    gridField(),
    bodySurfaceType_("bodySurface"),
    byName_(true),
    byPhysicalType_(false)
{
    nx_ = nx;
    ny_ = ny;
    nz_ = nz;

    xmin_ = x0;
    ymin_ = y0;
    zmin_ = z0;
    xmax_ = x1;
    ymax_ = y1;
    zmax_ = z1;

    generateNodes();
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::fieldSample::~fieldSample()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fieldSample::writeSamples(const word& fldName)
{
    const fileName path("./postProcessing");
    mkDir(path);

    word fname = "./postProcessing/" + fldName + ".dat";

    std::ofstream os(fname.c_str());
    if (!os)
    {
        FatalErrorInFunction
            << "Error writing sample field, file: \n"<<fname<<" cannot open."
            << exit(FatalError);
    }

    word content = "levelSet\n";

    std::vector<scalar> dvalues(nNodes(), 0);

    for (label i = 0; i < nx_ + 1; i++)
    {
        for (label j = 0; j < ny_ + 1; j++)
        {
            for (label k = 0; k < nzpt(); k++)
            {
                scalar val = gridField_(i, j, k);
                label id = ijk(i, j, k);

                dvalues[id] = val;
            }
        }
    }

    for (label ik = 0; ik < nNodes(); ik++)
    {
        content = content + helyxMap::doubleToString(dvalues[ik]) + "\n";
    }

    os  << content;
    os.close();
}


void Foam::fieldSample::sampleField
(
    const volScalarField& fld,
    const fvMesh* mesh,
    bool setLevel
)
{
    field_ = &fld;

    gridField_.setSize(nx_ + 1, ny_ + 1, nz_ + 1);

    inSolids_.setSize(nx_ + 1, ny_ + 1, nz_ + 1);

    for (label i = 0; i < nx_ + 1; i++)
    {
        for (label j = 0; j < ny_ + 1; j++)
        {
            for (label k = 0; k < nzpt(); k++)
            {
                point pt = gridNodes_(i, j, k);

                scalar xb = xbar(pt.x());
                scalar yb = ybar(pt.y());
                scalar zb = ybar(pt.z());

                std::vector<scalar> qvect(3);
                qvect[0] = xb;
                qvect[1] = yb;
                qvect[2] = zb;

                std::vector<label> nb(2);
                std::vector<scalar> dist(2);

                internalMap_.search(qvect, nb, dist);

                // nb stores the nearest cell index.
                // Search for the nearest body surface face center
                std::vector<label> nb1;
                std::vector<scalar> dist1;

                bodyMap_.search(qvect, nb1, dist1);

                // Check if a grid node is in solid
                label inSolid = 1;

                if (nb1.size() == 0)
                {
                    inSolid = false;    // cannot find boundary boundary match
                }
                else
                {
                    if (dist1[0] > dist[0])
                    {
                        inSolid = 0;    // nearer to fluid cell
                    }
                    else
                    {
                        // dist1[1] <= dist[0] can also be a fluid cell
                        // if cell_alpha > 1.0e-3
                        label bfaceid = nb1[0];

                        label r = bodyFaces_[bfaceid].first;
                        label jf = bodyFaces_[bfaceid].second;

                        point pf = mesh->Cf().boundaryField()[r][jf];

                        label ic = nb[0];
                        point pc = mesh->C()[ic];

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

                inSolids_(i, j, k) = inSolid;

                if (inSolid == 1)
                {
                    getSolidField(i, j, k, nb1, dist1, setLevel);
                }
                else
                {
                    getGridField(i, j, k, nb, dist);
                }
            }
        }
    }
}


void Foam::fieldSample::constructBodyKnn(const fvMesh* mesh)
{
    farray2d trainData;

    forAll(mesh->boundaryMesh(), r)
    {
        word type = mesh->boundaryMesh()[r].type();

        if (type == "empty" || type == "symmetryPlane")
        {
            continue;
        }

        word ptype = mesh->boundaryMesh()[r].physicalType();
        word bname = mesh->boundaryMesh()[r].name();

        if (!byName_)
        {
            if (!byPhysicalType_)
            {
                ptype = mesh->boundaryMesh()[r].type();
            }

            if (ptype != bodySurfaceType_)
            {
                continue;
            }
        }
        else
        {
            if (bname != bodySurfaceName_) continue;
        }

        forAll(mesh->boundaryMesh()[r], j)
        {
            point pt = mesh->Cf().boundaryField()[r][j];

            scalar xb = xbar(pt.x());
            scalar yb = ybar(pt.y());
            scalar zb = zbar(pt.z());

            std::vector<scalar> pts(3);
            pts[0] = xb;
            pts[1] = yb;
            pts[2] = zb;

            std::pair<label, label> faceid(r, j);

            bodyFaces_.append(faceid);

            trainData.push_back(pts);
        }
    }

    if (trainData.size() < 2)
    {
        FatalErrorInFunction
            << "Error constructing solid object search tree, "
            << "zero number of solid object boundary faces obtained." << nl
            << "Please check whether 'bodySurfaceType' name is set correctly."
            << exit(FatalError);
    }

    bodyMap_.build(trainData);
}


void Foam::fieldSample::constructFieldKnn(const fvMesh* mesh)
{
    farray2d trainData;

    forAll(mesh->C(), i)
    {
        point pt = mesh->C()[i];

        scalar xb = xbar(pt.x());
        scalar yb = ybar(pt.y());
        scalar zb = zbar(pt.z());

        std::vector<scalar> pts(3);
        pts[0] = xb;
        pts[1] = yb;
        pts[2] = zb;

        trainData.push_back(pts);
    }

    Info<< "Train data size: " << trainData.size() << endl;

    internalMap_.build(trainData);

    Info<< "Internal map build." << endl;
}


// ************************************************************************* //
