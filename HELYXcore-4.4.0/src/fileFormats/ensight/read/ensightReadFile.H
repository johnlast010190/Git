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
    (c) 2016 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::ensightReadFile

Description
    Ensight output with specialized read() for strings, integers and floats.
    Correctly handles binary read as well.

\*---------------------------------------------------------------------------*/

#ifndef ensightReadFile_H
#define ensightReadFile_H

#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/IOstreams/IOstreams/IOstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ensightReadFile Declaration
\*---------------------------------------------------------------------------*/

class ensightReadFile
:
    public IFstream
{
public:

    // Constructors

        //- Construct from pathname
        ensightReadFile
        (
            const fileName& pathname,
            IOstream::streamFormat format=IOstream::BINARY
        );

        //- Disallow default copy constructor
        ensightReadFile(const ensightReadFile&) = delete;


    //- Destructor
    ~ensightReadFile();


    // Output

        //- Inherit read from Istream
        using Istream::read;

        //- Binary read
        virtual Istream& read(char* buf, std::streamsize count);

        //- Read string as "%80s" or as binary
        Istream& read(string& value);

        //- Read integer as "%10d" or as binary
        Istream& read(label& value);

        //- Read float as "%12.5e" or as binary
        Istream& read(scalar& value);

        //- Read element keyword
        virtual Istream& readKeyword(string& key);

        //- Read "C Binary" for binary files (eg, geometry/measured)
        Istream& readBinaryHeader();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ensightReadFile&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
