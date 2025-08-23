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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2016-2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "MeshedSurfaceAllocator/MeshedSurfaceIOAllocator.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MeshedSurfaceIOAllocator::MeshedSurfaceIOAllocator
(
    const IOobject& ioPoints,
    const IOobject& ioFaces,
    const IOobject& ioZones
)
:
    points_(ioPoints),
    faces_(ioFaces),
    zones_(ioZones)
{}


Foam::MeshedSurfaceIOAllocator::MeshedSurfaceIOAllocator
(
    const IOobject& ioPoints,
    const pointField& points,
    const IOobject& ioFaces,
    const faceList& faces,
    const IOobject& ioZones,
    const surfZoneList& zones
)
:
    points_(ioPoints, points),
    faces_(ioFaces, faces),
    zones_(ioZones, zones)
{}


Foam::MeshedSurfaceIOAllocator::MeshedSurfaceIOAllocator
(
    const IOobject& ioPoints,
    pointField&& points,
    const IOobject& ioFaces,
    faceList&& faces,
    const IOobject& ioZones,
    surfZoneList&& zones
)
:
    points_(ioPoints, std::move(points)),
    faces_(ioFaces, std::move(faces)),
    zones_(ioZones, std::move(zones))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::MeshedSurfaceIOAllocator::~MeshedSurfaceIOAllocator()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::MeshedSurfaceIOAllocator::setInstance(const fileName& inst)
{
    points_.instance() = inst;
    faces_.instance()  = inst;
    zones_.instance()  = inst;
}


void Foam::MeshedSurfaceIOAllocator::setWriteOption(IOobject::writeOption w)
{
    points_.writeOpt() = w;
    faces_.writeOpt()  = w;
    zones_.writeOpt()  = w;
}


void Foam::MeshedSurfaceIOAllocator::clear()
{
    points_.clear();
    faces_.clear();
    zones_.clear();
}


void Foam::MeshedSurfaceIOAllocator::resetFaces
(
    List<face>&& faces,
    surfZoneList&& zones
)
{
    if (notNull(faces))
    {
        faces_.transfer(faces);
    }

    if (notNull(zones))
    {
        zones_.transfer(zones);
    }
}


void Foam::MeshedSurfaceIOAllocator::reset
(
    pointField&& points,
    faceList&& faces,
    surfZoneList&& zones
)
{
    // Take over new primitive data.
    // Optimized to avoid overwriting data at all
    if (notNull(points))
    {
        points_.transfer(points);
    }

    resetFaces(std::move(faces), std::move(zones));
}


void Foam::MeshedSurfaceIOAllocator::reset
(
    List<point>&& points,
    faceList&& faces,
    surfZoneList&& zones
)
{
    // Take over new primitive data.
    // Optimized to avoid overwriting data at all
    if (notNull(points))
    {
        points_.transfer(points);
    }

    resetFaces(std::move(faces), std::move(zones));
}


bool Foam::MeshedSurfaceIOAllocator::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    return
    (
        points_.writeObject(fmt, ver, cmp, write)
     && faces_.writeObject(fmt, ver, cmp, write)
     && zones_.writeObject(fmt, ver, cmp, write)
    );
}


// ************************************************************************* //
