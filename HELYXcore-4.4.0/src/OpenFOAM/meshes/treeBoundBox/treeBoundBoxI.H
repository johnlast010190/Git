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
    (c) 2017 OpenCFD Ltd.
    (c) 2010-2024 Engys Ltd.
    (c) 2011-2018 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "meshes/treeBoundBox/treeBoundBox.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::treeBoundBox::treeBoundBox()
:
    boundBox()
{}


inline Foam::treeBoundBox::treeBoundBox(const boundBox& bb)
:
    boundBox(bb)
{}


inline Foam::treeBoundBox::treeBoundBox(const point& pt)
:
    boundBox(pt)
{}


inline Foam::treeBoundBox::treeBoundBox(const point& min, const point& max)
:
    boundBox(min, max)
{}


inline Foam::treeBoundBox::treeBoundBox(Istream& is)
:
    boundBox(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::treeBoundBox::typDim() const
{
    return avgDim();
}


inline Foam::point Foam::treeBoundBox::corner(const direction octant) const
{
    return point
    (
        (octant & RIGHTHALF) ? max().x() : min().x(),
        (octant & TOPHALF)   ? max().y() : min().y(),
        (octant & FRONTHALF) ? max().z() : min().z()
    );
}


inline Foam::direction Foam::treeBoundBox::subOctant(const point& pt) const
{
    return subOctant(midpoint(), pt);
}


inline Foam::direction Foam::treeBoundBox::subOctant
(
    const point& mid,
    const point& pt
)
{
    direction octant = 0;

    if (pt.x() > mid.x())
    {
        octant |= treeBoundBox::RIGHTHALF;
    }

    if (pt.y() > mid.y())
    {
        octant |= treeBoundBox::TOPHALF;
    }

    if (pt.z() > mid.z())
    {
        octant |= treeBoundBox::FRONTHALF;
    }

    return octant;
}


inline Foam::direction Foam::treeBoundBox::subOctant
(
    const point& pt,
    bool& onEdge
) const
{
    return subOctant(midpoint(), pt, onEdge);
}


inline Foam::direction Foam::treeBoundBox::subOctant
(
    const point& mid,
    const point& pt,
    bool& onEdge
)
{
    direction octant = 0;
    onEdge = false;

    if (pt.x() > mid.x())
    {
        octant |= treeBoundBox::RIGHTHALF;
    }
    else if (pt.x() == mid.x())
    {
        onEdge = true;
    }

    if (pt.y() > mid.y())
    {
        octant |= treeBoundBox::TOPHALF;
    }
    else if (pt.y() == mid.y())
    {
        onEdge = true;
    }

    if (pt.z() > mid.z())
    {
        octant |= treeBoundBox::FRONTHALF;
    }
    else if (pt.z() == mid.z())
    {
        onEdge = true;
    }

    return octant;
}


inline Foam::direction Foam::treeBoundBox::subOctant
(
    const point& mid,
    const vector& dir,
    const point& pt,
    bool& onEdge
)
{
    direction octant = 0;
    onEdge = false;

    if (pt.x() > mid.x())
    {
        octant |= treeBoundBox::RIGHTHALF;
    }
    else if (pt.x() == mid.x())
    {
        onEdge = true;
        if (dir.x() > 0)
        {
            octant |= treeBoundBox::RIGHTHALF;
        }
    }

    if (pt.y() > mid.y())
    {
        octant |= treeBoundBox::TOPHALF;
    }
    else if (pt.y() == mid.y())
    {
        onEdge = true;
        if (dir.y() > 0)
        {
            octant |= treeBoundBox::TOPHALF;
        }
    }

    if (pt.z() > mid.z())
    {
        octant |= treeBoundBox::FRONTHALF;
    }
    else if (pt.z() == mid.z())
    {
        onEdge = true;
        if (dir.z() > 0)
        {
            octant |= treeBoundBox::FRONTHALF;
        }
    }

    return octant;
}


inline void Foam::treeBoundBox::searchOrder
(
    const point& pt,
    FixedList<direction,8>& octantOrder
) const
{
    vector dist = midpoint() - pt;

    direction octant = 0;

    if (dist.x() < 0)
    {
        octant |= treeBoundBox::RIGHTHALF;
        dist.x() *= -1;
    }

    if (dist.y() < 0)
    {
        octant |= treeBoundBox::TOPHALF;
        dist.y() *= -1;
    }

    if (dist.z() < 0)
    {
        octant |= treeBoundBox::FRONTHALF;
        dist.z() *= -1;
    }

    direction min = 0;
    direction mid = 0;
    direction max = 0;

    if (dist.x() < dist.y())
    {
        if (dist.y() < dist.z())
        {
            min = treeBoundBox::RIGHTHALF;
            mid = treeBoundBox::TOPHALF;
            max = treeBoundBox::FRONTHALF;
        }
        else if (dist.z() < dist.x())
        {
            min = treeBoundBox::FRONTHALF;
            mid = treeBoundBox::RIGHTHALF;
            max = treeBoundBox::TOPHALF;
        }
        else
        {
            min = treeBoundBox::RIGHTHALF;
            mid = treeBoundBox::FRONTHALF;
            max = treeBoundBox::TOPHALF;
        }
    }
    else
    {
        if (dist.z() < dist.y())
        {
            min = treeBoundBox::FRONTHALF;
            mid = treeBoundBox::TOPHALF;
            max = treeBoundBox::RIGHTHALF;
        }
        else if (dist.x() < dist.z())
        {
            min = treeBoundBox::TOPHALF;
            mid = treeBoundBox::RIGHTHALF;
            max = treeBoundBox::FRONTHALF;
        }
        else
        {
            min = treeBoundBox::TOPHALF;
            mid = treeBoundBox::FRONTHALF;
            max = treeBoundBox::RIGHTHALF;
        }
    }

    // Primary subOctant
    octantOrder[0] = octant;
    // subOctants joined to the primary by faces.
    octantOrder[1] = octant ^ min;
    octantOrder[2] = octant ^ mid;
    octantOrder[3] = octant ^ max;
    // subOctants joined to the primary by edges.
    octantOrder[4] = octantOrder[1] ^ mid;
    octantOrder[5] = octantOrder[1] ^ max;
    octantOrder[6] = octantOrder[2] ^ max;
    // subOctants joined to the primary by corner.
    octantOrder[7] = octantOrder[4] ^ max;
}


inline Foam::treeBoundBox Foam::treeBoundBox::extend
(
    const scalar s
) const
{
    // Numbers that don't approximate rational fractions with which to make the
    // box asymmetric. These are between zero and one.
    static const vector a((sqrt(5.0) + 1)/4, sqrt(2.0)/2, (sqrt(13.0) - 1)/4);
    static const vector b(a.y(), a.z(), a.x());

    treeBoundBox bb(*this);

    const scalar delta = s*Foam::mag(bb.span());
    bb.min() -= Foam::max(delta*a, vector::uniform(ROOTVSMALL));
    bb.max() += Foam::max(delta*b, vector::uniform(ROOTVSMALL));

#if defined(HELYX_SP)
    // extend in single precision if needed
    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        if (bb.min()[dir] >= bb.max()[dir])
        {
            bb.min()[dir] = nextafterf(bb.min()[dir],-INFINITY);
            bb.max()[dir] = nextafterf(bb.max()[dir],INFINITY);
        }
    }
#endif

    return bb;
}

// * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * * //

inline bool Foam::operator==(const treeBoundBox& a, const treeBoundBox& b)
{
    return static_cast<const boundBox&>(a) == static_cast<const boundBox&>(b);
}


inline bool Foam::operator!=(const treeBoundBox& a, const treeBoundBox& b)
{
    return !(a == b);
}

// ************************************************************************* //
