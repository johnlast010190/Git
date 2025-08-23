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
    (c) 2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "db/IOstreams/Fstreams/masterOFstream.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "include/OSspecific.H"
#include "db/IOstreams/Pstreams/PstreamBuffers.H"
#include "global/fileOperations/masterUncollatedFileOperation/masterUncollatedFileOperation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::masterOFstream::checkWrite
(
    const fileName& fName,
    const string& str
)
{
    mkDir(fName.path());

    OFstream os
    (
        fName,
        IOstream::BINARY,   //format(),
        version(),
        compression_,
        append_
    );
    if (!os.good())
    {
        FatalIOErrorInFunction(os)
            << "Could not open file " << fName
            << exit(FatalIOError);
    }

    os.writeQuoted(str, false);
    if (!os.good())
    {
        FatalIOErrorInFunction(os)
            << "Failed writing to " << fName
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::masterOFstream::masterOFstream
(
    const fileName& pathName,
    streamFormat format,
    versionNumber version,
    compressionType compression,
    const bool append,
    const bool write
)
:
    OStringStream(format, version),
    pathName_(pathName),
    compression_(compression),
    append_(append),
    write_(write)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::masterOFstream::~masterOFstream()
{
    if (Pstream::parRun())
    {
        List<fileName> filePaths(Pstream::nProcs());
        filePaths[Pstream::myProcNo()] = pathName_;
        Pstream::gatherList(filePaths);

        bool uniform =
            fileOperations::masterUncollatedFileOperation::uniformFile
            (
                filePaths
            );

        Pstream::scatter(uniform);

        if (uniform)
        {
            if (Pstream::master() && write_)
            {
                checkWrite(pathName_, str());
            }
            return;
        }
        boolList write(Pstream::nProcs());
        write[Pstream::myProcNo()] = write_;
        Pstream::gatherList(write);

        // Different files
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send my buffer to master
        if (!Pstream::master())
        {
            UOPstream os(Pstream::masterNo(), pBufs);
            string s(this->str());
            os.write(&s[0], s.size());
        }

        labelList recvSizes;
        pBufs.finishedSends(recvSizes);

        if (Pstream::master())
        {
            // Write my own data
            {
                if (write[Pstream::myProcNo()])
                {
                    checkWrite(filePaths[Pstream::myProcNo()], str());
                }
            }

            for (label proci = 1; proci < Pstream::nProcs(); proci++)
            {
                UIPstream is(proci, pBufs);
                List<char> buf(recvSizes[proci]);

                is.read(buf.begin(), buf.size());

                if (write[proci])
                {
                    checkWrite
                    (
                        filePaths[proci],
                        string(buf.begin(), buf.size())
                    );
                }
            }
        }
    }
    else
    {
        checkWrite(pathName_, str());
    }
}


// ************************************************************************* //
