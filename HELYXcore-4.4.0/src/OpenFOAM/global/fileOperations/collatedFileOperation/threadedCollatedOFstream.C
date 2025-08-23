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
    (c) 2017-2018 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "global/fileOperations/collatedFileOperation/threadedCollatedOFstream.H"
#include "db/IOobjects/decomposedBlockData/decomposedBlockData.H"
#include "global/fileOperations/collatedFileOperation/OFstreamCollator.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::threadedCollatedOFstream::threadedCollatedOFstream
(
    OFstreamCollator& writer,
    const fileName& pathName,
    streamFormat format,
    versionNumber version,
    compressionType compression,
    const bool useThread
)
:
    OStringStream(format, version),
    writer_(writer),
    pathName_(pathName),
    compression_(compression),
    useThread_(useThread)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::threadedCollatedOFstream::~threadedCollatedOFstream()
{
    writer_.write
    (
        decomposedBlockData::typeName,
        pathName_,
        str(),
        IOstream::BINARY,
        version(),
        compression_,
        false,                  // append
        useThread_
    );
}


// ************************************************************************* //
