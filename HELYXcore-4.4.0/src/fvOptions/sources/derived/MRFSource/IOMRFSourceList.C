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
    (c) 2012-2014 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "IOMRFSourceList.H"
#include "fvMesh/fvMesh.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::IOobject Foam::fv::IOMRFSourceList::createIOobject
(
    const objectRegistry& obr
) const
{
    IOobject io
    (
        "fvOptions",
        obr.time().system(),
        obr,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.typeHeaderOk<IOdictionary>(true))
    {
        Info<< "Creating MRF sources list from " << io.name() << nl << endl;

        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        Info<< "No MRF models present" << nl << endl;

        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::IOMRFSourceList::IOMRFSourceList
(
    const objectRegistry& obr
)
:
    IOdictionary(createIOobject(obr)),
    MRFSourceList(obr, *this)
{}


bool Foam::fv::IOMRFSourceList::read()
{
    if (regIOobject::read())
    {
        MRFSourceList::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
