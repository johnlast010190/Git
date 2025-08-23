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
    (c) 2020-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "meshes/smartBoundBox/smartBoundBox.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::smartBoundBox::smartBoundBox()
:
    boundBox(),
    boxType_(0),
    ampCoef_(1.5),
    nx_(100),
    ny_(100),
    nz_(1)
{}


Foam::smartBoundBox::smartBoundBox
(
    const point& pmin,
    const point& pmax,
    const word type,
    scalar ampCoef
)
:
    boundBox(pmin, pmax),
    nx_(100),
    ny_(100),
    nz_(1)
{
    if (type == "boundBox")
    {
        boxType_ = 0;
        ampCoef_ = 1;
    }
    else
    {
        boxType_ = 1;
        ampCoef_ = ampCoef;
    }
}


Foam::smartBoundBox::smartBoundBox
(
    const fvMesh* mesh,
    const word type,
    scalar ampCoef,
    label cellNum
)
:
    boundBox(mesh->C(), false),
    ampCoef_(ampCoef),
    nx_(100),
    ny_(100),
    nz_(1)
{
    if (type == "simple")
    {
        boxType_ = 0;
    }
    else
    {
        boxType_ = 1;
    }

    createBox(mesh, cellNum);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::smartBoundBox::getCellList(const fvMesh* mesh, label nCells)
{
    const point& pmin = boundBox::min();
    const point& pmax = boundBox::max();

    Info<< "\nsmartBoundBox:" << "\tmin: " << pmin << "\tmax: " << pmax << endl;
    Info<< "nCells: " << nCells << nl << endl;

    scalar deltaX = pmax.x() - pmin.x();
    scalar deltaY = pmax.y() - pmin.y();
    scalar deltaZ = pmax.z() - pmin.z();

    Info<< "deltaXYZ: " << deltaX << " " << deltaY << " " << deltaZ << endl;

    scalar dmax = Foam::max(deltaX, deltaY);
    dmax = Foam::max(dmax, deltaZ);

    scalar ax = deltaX/(dmax + SMALL);
    scalar ay = deltaY/(dmax + SMALL);
    scalar az = deltaZ/(dmax + SMALL);

    nx_ = 1;
    ny_ = 1;
    nz_ = 1;

    scalar ann = sqrt(scalar(nCells));

    if (ax*ann < 1.0 && ay*az > 0.0)
    {
        // 2d in x direction
        nz_ = label(sqrt(az*nCells/ay) + 0.5);
        ny_ = label((ay/az)*nz_ + 0.5);
        nx_ = 1;
    }
    else if (ay*ann < 1.0 && ax*az > 0.0)
    {
        // 2d in y direction
        nz_ = label(sqrt(az*nCells/ax) + 0.5);
        nx_ = label((ax/az)*nz_ + 0.5);
        ny_ = 1;
    }
    else if (az*ann < 1.0 && ax*ay > 0.0)
    {
        // 2d in z direction
        ny_ = label(sqrt(ay*nCells/ax) + 0.5);
        nx_ = label((ax/ay)*ny_ + 0.5);
        nz_ = 1;
    }
    else if (ax*ay*az > 0.0)
    {
        scalar coef = az*az/(ax*ay);
        nz_ = label(Foam::pow(nCells*coef, 1.0/3.0) + 0.5);
        nx_ = label((ax/az)*nz_ + 0.5);
        ny_ = label((ay/az)*nz_ + 0.5);
    }

    Info<< "nx, ny, nz: " << nx_ << " " << ny_ << " " << nz_ << endl;

    numCells_.setSize(nx_ + 1, ny_ + 1, nz_ + 1);
    numCells_ = 0;

    deltaX_ = deltaX/scalar(nx_) + SMALL;
    deltaY_ = deltaY/scalar(ny_) + SMALL;
    deltaZ_ = deltaZ/scalar(nz_) + SMALL;

    Info<< "deltas: " << deltaX_ << " " << deltaY_ << " "
        << deltaZ_ << nl << endl;

    ampCoef_ = Foam::max(ampCoef_, 1);
    ampCoef_ = Foam::min(ampCoef_, 1.5);

    // Determine the number of cells in each grid
    forAll(mesh->C(), i)
    {
        const point& pt = mesh->C()[i];

        if (!contains(pt)) continue;

        label ix = label(Foam::min(labelMax, (pt.x() - pmin.x())/deltaX_));
        label iy = label(Foam::min(labelMax, (pt.y() - pmin.y())/deltaY_));
        label iz = label(Foam::min(labelMax, (pt.z() - pmin.z())/deltaZ_));

        ix = Foam::min(nx_ - 1, ix);
        iy = Foam::min(ny_ - 1, iy);
        iz = Foam::min(nz_ - 1, iz);

        numCells_.setVal(ix, iy, iz, true);
    }

    // Second level
    const vectorField verts = mesh->points();

    forAll(mesh->cellPoints(), ic)
    {
        const labelUList& cverts = mesh->cellPoints()[ic];

        scalar xminc =  GREAT;
        scalar xmaxc = -GREAT;
        scalar yminc =  GREAT;
        scalar ymaxc = -GREAT;
        scalar zminc =  GREAT;
        scalar zmaxc = -GREAT;

        forAll(cverts, j)
        {
            label jv = cverts[j];

            const point& pv = verts[jv];

            xminc = Foam::min(xminc, pv.x());
            yminc = Foam::min(yminc, pv.y());
            zminc = Foam::min(zminc, pv.z());

            xmaxc = Foam::max(xmaxc, pv.x());
            ymaxc = Foam::max(ymaxc, pv.y());
            zmaxc = Foam::max(zmaxc, pv.z());
        }

        xmaxc = xminc + ampCoef_*(xmaxc - xminc);
        ymaxc = yminc + ampCoef_*(ymaxc - yminc);
        zmaxc = zminc + ampCoef_*(zmaxc - zminc);

        label i0 =
            label
            (
                Foam::max(0, Foam::min(labelMax, (xminc - pmin.x())/deltaX_))
            );
        label i1 =
            label
            (
                Foam::min(nx_, Foam::max(labelMin, (xmaxc - pmin.x())/deltaX_))
            );

        label j0 =
            label
            (
                Foam::max(0, Foam::min(labelMax, (yminc - pmin.y())/deltaY_))
            );
        label j1 =
            label
            (
                Foam::min(ny_, Foam::max(labelMin, (ymaxc - pmin.y())/deltaY_))
            );

        label k0 =
            label
            (
                Foam::max(0, Foam::min(labelMax, (zminc - pmin.z())/deltaZ_))
            );
        label k1 =
            label
            (
                Foam::min(nz_, Foam::max(labelMin, (zmaxc - pmin.z())/deltaZ_))
            );

        for (label i = i0; i < i1; i++)
        {
            for (label j = j0; j < j1; j++)
            {
                for (label k = k0; k < k1; k++)
                {
                    if (!numCells_.getVal(i, j, k))
                    {
                        numCells_.setVal(i, j, k, true);
                    }
                }
            }
        }
    }
}


void Foam::smartBoundBox::createBox(const fvMesh* mesh, label cellNum)
{
    if (boxType_ == 1)
    {
        getCellList(mesh, cellNum);
    }
}


// ************************************************************************* //
