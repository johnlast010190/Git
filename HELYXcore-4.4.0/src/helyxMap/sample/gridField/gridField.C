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

#include "sample/gridField/gridField.H"
#include "helyxMap.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::gridField::rhoRef_;

Foam::vector Foam::gridField::Uref_;



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gridField::gridField()
:
    containSolidObject_(true),
    nx_(20),
    ny_(20),
    nz_(1),
    dx_(0),
    dy_(0),
    dz_(0),
    interp_(false),
    interp2_(false)
{}


Foam::gridField::gridField
(
    const boundBox& bbox,
    label nx,
    label ny,
    label nz,
    bool containObj
)
:
    containSolidObject_(containObj),
    interp_(false)
{
    getBbox(bbox);

    setDivision(nx, ny, nz);

    setDelta();
}


Foam::gridField::gridField
(
    const word& bbox,
    label nx,
    label ny,
    label nz,
    bool containObj
)
:
    containSolidObject_(containObj),
    interp_(false)
{
    getBbox(bbox);

    setDivision(nx, ny, nz);

    setDelta();
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::gridField::~gridField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::gridField::inList(const word& key, const wordList& vals)
{
    label sz = vals.size();

    for (label i = 0; i < sz; i++)
    {
        if (vals[i] == key)
        {
            return i;
        }
    }

    return -1;
}


void Foam::gridField::setDivision(label nx, label ny, label nz)
{
    // Called after set field names
    nx_ = nx;
    ny_ = ny;
    nz_ = nz;
}


void Foam::gridField::setDelta()
{
    dx_ = (xmax() - xmin())/scalar(nx_);
    dy_ = (ymax() - ymin())/scalar(ny_);
    dz_ = (zmax() - zmin())/scalar(nz_);
}


void Foam::gridField::setSampleNames(const word& scalars, const word& vectors)
{
    // Names are separated by comma
    helyxMap::stringTokenize(scalarSamples_, scalars, ",");

    if (vectors != "none")
    {
        helyxMap::stringTokenize(vectorSamples_, vectors, ",");
    }
}


void Foam::gridField::generateNodes()
{
    // Generate node coordinates and store in gridNodes
    dx_= (xmax() - xmin())/scalar(nx_);
    dy_= (ymax() - ymin())/scalar(ny_);
    dz_= (zmax() - zmin())/scalar(nzpt());

    if (dx_ < SMALL || dy_ < SMALL || dz_ < SMALL)
    {
        FatalErrorInFunction
            << "One or more of the dx,dy,dz are too small. "
            << "Please check whether bounding box is set properly." << nl
            << "dx = " << dx_ << ", dy = " << dy_ << ", dz = " << dz_
            << exit(FatalError);
    }

    if (nz_ > 1)
    {
        gridNodes_.setSize(nx_ + 1, ny_ + 1, nz_ + 1);
    }
    else
    {
        gridNodes_.setSize(nx_ + 1, ny_ + 1, 1);
    }

    for (label i = 0; i < nx_ + 1; ++i)
    {
        scalar xi = xmin() + i*dx_;

        for (label j = 0; j < ny_ + 1; ++j)
        {
            scalar yj = ymin() + j*dy_;

            for (label k = 0; k < nzpt(); ++k)
            {
                scalar zk = zmin() + k*dz_;

                point node(xi, yj, zk);
                gridNodes_(i, j, k) = node;
            }
        }
    }
}


void Foam::gridField::getBbox(const boundBox& bbox)
{
    xmin_ = bbox.min().x();
    ymin_ = bbox.min().y();
    zmin_ = bbox.min().z();
    xmax_ = bbox.max().x();
    ymax_ = bbox.max().y();
    zmax_ = bbox.max().z();
}


void Foam::gridField::getBbox(const word& val)
{
    wordList vtmp;

    helyxMap::stringTokenize(vtmp, val, ",");

    if (vtmp.size() != 6)
    {
        FatalErrorInFunction
            << "Incorrect data: " << val
            << exit(FatalError);
    }

    xmin_ = atof(vtmp[0].c_str());
    ymin_ = atof(vtmp[1].c_str());
    zmin_ = atof(vtmp[2].c_str());
    xmax_ = atof(vtmp[3].c_str());
    ymax_ = atof(vtmp[4].c_str());
    zmax_ = atof(vtmp[5].c_str());
}


void Foam::gridField::setSampleNames
(
    const HashSet<word>& scalars,
    const HashSet<word>& vectors
)
{
    forAllConstIter(HashSet<word>, scalars, iter)
    {
        scalarSamples_.append(iter());
    }

    forAllConstIter(HashSet<word>, vectors, iter)
    {
        vectorSamples_.append(iter());
    }

    setFieldSize();
}


void Foam::gridField::setFieldSize()
{
    if (gridNodes_.size() == 0)
    {
        gridNodes_.setSize(nx_ + 1, ny_ + 1, nzpt());
    }

    inSolids_.setSize(nx_ + 1, ny_ + 1, nzpt());

    if (scalarSamples_.size() >= 1)
    {
        gridFields_.setSize(nx_ + 1, ny_ + 1, nzpt());
    }

    if (vectorSamples_.size() >= 1)
    {
        gridVectors_.setSize(nx_ + 1, ny_ + 1, nzpt());
    }

    dists_.setSize(nx_ + 1, ny_ + 1, nzpt());
}


// ************************************************************************* //
