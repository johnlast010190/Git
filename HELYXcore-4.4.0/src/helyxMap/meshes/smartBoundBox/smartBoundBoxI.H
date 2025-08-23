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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::smartBoundBox::pointInBox(const point& pt) const
{
    bool inBoundBox = contains(pt);

    if (boxType_ == 0)
    {
        return inBoundBox;
    }

    if (!inBoundBox)
    {
        return false;
    }

    const point& pmin = boundBox::min();

    label ix = label((pt.x() - pmin.x())/deltaX_);
    label iy = label((pt.y() - pmin.y())/deltaY_);
    label iz = label((pt.z() - pmin.z())/deltaZ_);

    if (ix > nx_ - 1 || iy > ny_ - 1 || iz > nz_ - 1)
    {
        return false;
    }

    return numCells_(ix, iy, iz);
}


// ************************************************************************* //
