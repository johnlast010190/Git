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
    (c) 2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fvMeshDistribute/IOmapDistributePolyMesh.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
defineTypeNameAndDebug(IOpolyDistributionMap, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IOpolyDistributionMap::IOpolyDistributionMap(const IOobject& io)
:
    regIOobject(io)
{
    // Warn for MUST_READ_IF_MODIFIED
    warnNoRereading<IOpolyDistributionMap>();

    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }
}


Foam::IOpolyDistributionMap::IOpolyDistributionMap
(
    const IOobject& io,
    const polyDistributionMap& map
)
:
    regIOobject(io)
{
    // Warn for MUST_READ_IF_MODIFIED
    warnNoRereading<IOpolyDistributionMap>();

    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }
    else
    {
        polyDistributionMap::operator=(map);
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::IOpolyDistributionMap::~IOpolyDistributionMap()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::IOpolyDistributionMap::readData(Istream& is)
{
    return (is >> *this).good();
}


bool Foam::IOpolyDistributionMap::writeData(Ostream& os) const
{
    return (os << *this).good();
}


// ************************************************************************* //
