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

Class
    Foam::masterOFstream

Description
    Master-only drop-in replacement for OFstream.

SourceFiles
    masterOFstream.C

\*---------------------------------------------------------------------------*/

#ifndef masterOFstream_H
#define masterOFstream_H

#include "db/IOstreams/StringStreams/OStringStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class masterOFstream Declaration
\*---------------------------------------------------------------------------*/

class masterOFstream
:
    public OStringStream
{
    // Private data

        const fileName pathName_;

        const IOstream::compressionType compression_;

        const bool append_;

        //- Should file be written
        const bool write_;


    // Private Member Functions

        //- Open file with checking
        void checkWrite(const fileName& fName, const string& str);


public:

    // Constructors

        //- Construct and set stream status
        masterOFstream
        (
            const fileName& pathname,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            compressionType compression=UNCOMPRESSED,
            const bool append = false,
            const bool write = true
        );


    //- Destructor
    ~masterOFstream();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
