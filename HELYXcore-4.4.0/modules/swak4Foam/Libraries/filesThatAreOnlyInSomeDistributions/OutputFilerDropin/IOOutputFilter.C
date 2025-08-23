/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : dev
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
    (c) 2011-2012 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "IOOutputFilter.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class OutputFilter>
Foam::IOOutputFilter<OutputFilter>::IOOutputFilter
(
    const word& outputFilterName,
    const IOobject& ioDict,
    const bool readFromFiles
)
:
    IOdictionary(ioDict),
    OutputFilter(outputFilterName, ioDict.db(), *this, readFromFiles)
{}


template<class OutputFilter>
Foam::IOOutputFilter<OutputFilter>::IOOutputFilter
(
    const word& outputFilterName,
    const objectRegistry& obr,
    const word& dictName,
    const IOobject::readOption rOpt,
    const bool readFromFiles
)
:
    IOdictionary
    (
        IOobject
        (
            dictName,
            obr.time().system(),
            obr,
            rOpt,
            IOobject::NO_WRITE
        )
    ),
    OutputFilter(outputFilterName, obr, *this, readFromFiles)
{}


template<class OutputFilter>
Foam::IOOutputFilter<OutputFilter>::IOOutputFilter
(
    const word& outputFilterName,
    const objectRegistry& obr,
    const fileName& dictName,
    const IOobject::readOption rOpt,
    const bool readFromFiles
)
:
    IOdictionary
    (
        IOobject
        (
            dictName,
            obr,
            rOpt,
            IOobject::NO_WRITE
        )
    ),
    OutputFilter(outputFilterName, obr, *this, readFromFiles)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class OutputFilter>
Foam::IOOutputFilter<OutputFilter>::~IOOutputFilter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class OutputFilter>
bool Foam::IOOutputFilter<OutputFilter>::read()
{
    if (regIOobject::read())
    {
        OutputFilter::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}


template<class OutputFilter>
bool Foam::IOOutputFilter<OutputFilter>::write()
{
    OutputFilter::write();
    return true;
}


template<class OutputFilter>
void Foam::IOOutputFilter<OutputFilter>::movePoints(const polyMesh& mesh)
{
    read();
    OutputFilter::movePoints(mesh);
}


template<class OutputFilter>
void Foam::IOOutputFilter<OutputFilter>::topoChange
(
    const polyTopoChangeMap& map
)
{
    read();
    OutputFilter::topoChange(map);
}


template<class OutputFilter>
void Foam::IOOutputFilter<OutputFilter>::mapMesh
(
    const polyMeshMap& map
)
{
    read();
    OutputFilter::mapMesh(map);
}


// ************************************************************************* //
