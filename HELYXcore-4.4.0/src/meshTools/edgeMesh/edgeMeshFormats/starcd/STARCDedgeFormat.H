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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::fileFormats::STARCDedgeFormat

Description
    Read/write the lines from pro-STAR vrt/cel files.

Note
    Uses the extension \a .inp (input) to denote the format.

See also
    Foam::meshReaders::STARCD

SourceFiles
    STARCDedgeFormat.C

\*---------------------------------------------------------------------------*/

#ifndef STARCDedgeFormat_H
#define STARCDedgeFormat_H

#include "edgeMesh/edgeMesh.H"
#include "starcd/STARCDCore.H"

#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/IOstreams/IOstreams/Ostream.H"
#include "db/IOstreams/Fstreams/OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                      Class STARCDedgeFormat Declaration
\*---------------------------------------------------------------------------*/

class STARCDedgeFormat
:
    public edgeMesh,
    protected STARCDCore
{
    // Private Member Functions

        static inline void writeLines(Ostream&, const edgeList&);


protected:

    // Protected Member Functions

    static void writeCase
    (
        Ostream&,
        const pointField&,
        const label nEdges
    );


public:

    // Constructors

        //- Construct from file name
        STARCDedgeFormat(const fileName&);

        //- Disallow default bitwise copy construction
        STARCDedgeFormat(const STARCDedgeFormat&) = delete;


    // Selectors

        //- Read file and return edgeMesh
        static autoPtr<edgeMesh> New(const fileName& name)
        {
            return autoPtr<edgeMesh>
            (
                new STARCDedgeFormat(name)
            );
        }


    //- Destructor
    virtual ~STARCDedgeFormat()
    {}


    // Member Functions

        //- Write edge mesh
        static void write(const fileName&, const edgeMesh&);

        //- Read from file
        virtual bool read(const fileName&);

        //- Write object
        virtual void write(const fileName& name) const
        {
            write(name, *this);
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const STARCDedgeFormat&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
