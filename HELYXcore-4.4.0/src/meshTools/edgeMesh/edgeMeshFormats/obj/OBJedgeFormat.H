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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::fileFormats::OBJedgeFormat

Description
    Provide a means of reading/writing Alias/Wavefront OBJ format.

    Does not handle negative vertex indices.

SourceFiles
    OBJedgeFormat.C

\*---------------------------------------------------------------------------*/

#ifndef OBJedgeFormat_H
#define OBJedgeFormat_H

#include "edgeMesh/edgeMesh.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/IOstreams/IOstreams/Ostream.H"
#include "db/IOstreams/Fstreams/OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                        Class OBJedgeFormat Declaration
\*---------------------------------------------------------------------------*/

class OBJedgeFormat
:
    public edgeMesh
{
    // Private Member Functions

        void readVertices
        (
            const string& line,
            string::size_type& endNum,
            DynamicList<label>& dynVertices
        );


public:

    // Constructors

        //- Construct from file name
        OBJedgeFormat(const fileName&);

        //- Disallow default bitwise copy construction
        OBJedgeFormat(const OBJedgeFormat&) = delete;


    // Selectors

        //- Read file and return surface
        static autoPtr<edgeMesh> New(const fileName& name)
        {
            return autoPtr<edgeMesh>
            (
                new OBJedgeFormat(name)
            );
        }


    //- Destructor
    virtual ~OBJedgeFormat()
    {}


    // Member Functions

        //- Write surface mesh components by proxy
        static void write(const fileName&, const edgeMesh&);

        //- Read from file
        virtual bool read(const fileName&);

        //- Write object file
        virtual void write(const fileName& name) const
        {
            write(name, *this);
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const OBJedgeFormat&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
