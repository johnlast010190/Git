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
    (c) 2011-2014 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::OSstream

Description
    Generic output stream.

SourceFiles
    OSstreamI.H
    OSstream.C

\*---------------------------------------------------------------------------*/

#ifndef OSstream_H
#define OSstream_H

#include "db/IOstreams/IOstreams/Ostream.H"

#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class OSstream Declaration
\*---------------------------------------------------------------------------*/

class OSstream
:
    public Ostream
{
    // Private data

        fileName name_;
        ostream& os_;


public:

    // Constructors

        //- Set stream status
        OSstream
        (
            ostream& os,
            const string& name,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            compressionType compression=UNCOMPRESSED
        );


    // Member functions

        // Enquiry

            //- Return the name of the stream
            //  Useful for Fstream to return the filename
            virtual const fileName& name() const
            {
                return name_;
            }

            //- Return non-const access to the name of the stream
            //  Useful to alter the stream name
            virtual fileName& name()
            {
                return name_;
            }

            //- Return flags of output stream
            virtual ios_base::fmtflags flags() const;


        // Write functions

            //- Write next token to stream
            virtual Ostream& write(const token&);

            //- Write character
            virtual Ostream& write(const char);

            //- Write character string
            virtual Ostream& write(const char*);

            //- Write word
            virtual Ostream& write(const word&);

            //- Write string
            //  In the rare case that the string contains a final trailing
            //  backslash, it will be dropped to the appearance of an escaped
            //  double-quote.
            virtual Ostream& write(const string&);

            //- Write std::string surrounded by quotes.
            //  Optional write without quotes.
            virtual Ostream& writeQuoted
            (
                const std::string&,
                const bool quoted=true
            );

            //- Write int32_t
            virtual Ostream& write(const int32_t);

            //- Write int64_t
            virtual Ostream& write(const int64_t);

            //- Write floatScalar
            virtual Ostream& write(const floatScalar);

            //- Write doubleScalar
            virtual Ostream& write(const doubleScalar);

            //- Write binary block
            virtual Ostream& write(const char*, std::streamsize);

            //- Add indentation characters
            virtual void indent();


        // Stream state functions

            //- Set flags of output stream
            virtual ios_base::fmtflags flags(const ios_base::fmtflags flags);

            //- Flush stream
            virtual void flush();

            //- Add newline and flush stream
            virtual void endl();

            //- Get width of output field
            virtual int width() const;

            //- Set width of output field (and return old width)
            virtual int width(const int);

            //- Get precision of output field
            virtual int precision() const;

            //- Set precision of output field (and return old precision)
            virtual int precision(const int);


        // STL stream

            //- Access to underlying std::ostream
            virtual ostream& stdStream()
            {
                return os_;
            }

            //- Const access to underlying std::ostream
            virtual const ostream& stdStream() const
            {
                return os_;
            }


        // Print

            //- Print description of IOstream to Ostream
            virtual void print(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const OSstream&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "db/IOstreams/Sstreams/OSstreamI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
