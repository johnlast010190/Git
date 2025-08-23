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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::point Foam::gridField::gridPoint(label i, label j, label k)
{
    point pt(xmin() + dx_*i, ymin() + dy_*j, zmin() + dz_*k);

    return pt;
}


inline Foam::scalar Foam::gridField::p(label i, label j, label k, label ip)
{
    return gridFields_(i, j, k)[ip];
}


inline Foam::point Foam::gridField::U(label i, label j, label k, label iv)
{
    return gridVectors_(i, j, k)[iv];
}


inline void Foam::gridField::getGridId
(
    label& i,
    label& j,
    label& k,
    scalar x,
    scalar y,
    scalar z
)
{
    i = label((x - xmin())/dx_);
    j = label((y - ymin())/dy_);
    k = 0;

    if (dimension() == 3)
    {
        k = label((z - zmin())/dz_);
    }
}


inline Foam::label Foam::gridField::scalarId(const word& scname)
{
    for (label isc = 0; isc < scalarSamples_.size(); isc++)
    {
        if (scalarSamples_[isc] == scname)
        {
            return isc;
        }
    }

    return -1;
}


inline Foam::point Foam::gridField::getVector
(
    word vecName,
    scalar x,
    scalar y,
    scalar z
)
{
    label iv = vectorId(vecName);

    if (iv < 0)
    {
        FatalErrorInFunction
            << "Error getting vector from grid, " << vecName << " not found.\n"
            << exit(FatalError);
    }

    label i, j, k;
    getGridId(i, j, k, x, y, z);

    point v1 = gridVectors_(i, j, k)[iv];

    if (dimension() == 3 || !interp2_)
    {
        return v1;
    }

    label i1 = i + 1;
    if (i1 > nx_)
    {
        i1 = nx_;
    }

    label j1 = j + 1;
    if (j1 > ny_)
    {
        j1 = ny_;
    }

    point v2 = gridVectors_(i1, j, k)[iv];
    point v3 = gridVectors_(i1, j1, k)[iv];
    point v4 = gridVectors_(i, j1, k)[iv];

    point p0(x, y, z);
    point p1(xmin() + dx_*i, ymin() + dy_*j, z);
    point p2(xmin() + dx_*i1, ymin() + dy_*j, z);
    point p3(xmin() + dx_*i1, ymin() + dy_*j1, z);
    point p4(xmin() + dx_*i, ymin() + dy_*j1, z);

    point p01 = p0 - p1;
    point p02 = p0 - p2;
    point p03 = p0 - p3;
    point p04 = p0 - p4;

    scalar s01 = mag(p01) + 1.0e-20;
    scalar s02 = mag(p02) + 1.0e-20;
    scalar s03 = mag(p03) + 1.0e-20;
    scalar s04 = mag(p04) + 1.0e-20;

    point val = v1/s01 + v2/s02 + v3/s03 + v4/s04;

    scalar ds = 1/s01 + 1/s02 + 1/s03 + 1/s04;

    return val/ds;
}


inline Foam::scalar Foam::gridField::getScalar
(
    word scName,
    scalar x,
    scalar y,
    scalar z
)
{
    label isc = scalarId(scName);

    if (isc < 0)
    {
        FatalErrorInFunction
            << "Error getting scalar from grid, " << scName << " not found.\n"
            << exit(FatalError);
    }

    label i, j, k;
    getGridId(i, j, k, x, y, z);

    scalar s1 = gridFields_(i, j, k)[isc];

    if (dimension() == 3 || !interp2_)
    {
        return s1;
    }

    label i1 = i + 1;
    if (i1 > nx_)
    {
        i1 = nx_;
    }

    label j1 = j + 1;
    if (j1 > ny_)
    {
        j1 = ny_;
    }

    scalar s2 = gridFields_(i1, j, k)[isc];
    scalar s3 = gridFields_(i1, j1, k)[isc];
    scalar s4 = gridFields_(i, j1, k)[isc];

    point p0(x, y, z);
    point p1(xmin() + dx_*i, ymin() + dy_*j, z);
    point p2(xmin() + dx_*i1, ymin() + dy_*j, z);
    point p3(xmin() + dx_*i1, ymin() + dy_*j1, z);
    point p4(xmin() + dx_*i, ymin() + dy_*j1, z);

    point p01 = p0 - p1;
    point p02 = p0 - p2;
    point p03 = p0 - p3;
    point p04 = p0 - p4;

    scalar s01 = mag(p01) + 1.0e-20;
    scalar s02 = mag(p02) + 1.0e-20;
    scalar s03 = mag(p03) + 1.0e-20;
    scalar s04 = mag(p04) + 1.0e-20;

    scalar val = s1/s01 + s2/s02 + s3/s03 + s4/s04;

    scalar ds = 1/s01 + 1/s02 + 1/s03 + 1/s04;

    return val/ds;
}


inline Foam::label Foam::gridField::vectorId(const word& vecname)
{
    for (label iv = 0; iv < vectorSamples_.size(); iv++)
    {
        if (vectorSamples_[iv] == vecname)
        {
            return iv;
        }
    }

    return -1;
}


inline bool Foam::gridField::inBlock(scalar x, scalar y, scalar z)
{
    if
    (
        x >= xmin() && x <= xmax()
     && y >= ymin() && y <= ymax()
     && z >= zmin() && z <= zmax()
    )
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
