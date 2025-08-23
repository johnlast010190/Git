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
    (c) 2017-2024 Engys Ltd.

Class
    Foam::fileFormats::edgeMeshFormat

Description
    Provide a means of reading/writing the single-file OpenFOAM edge format.

Note
   This class provides more methods than the regular edge format interface.

SourceFiles
    edgeMeshFormat.C

\*---------------------------------------------------------------------------*/

#ifndef edgeMeshFormat_H
#define edgeMeshFormat_H

#include "edgeMesh/edgeMesh.H"
#include "db/IOstreams/IOstreams/Ostream.H"
#include "db/IOstreams/Fstreams/OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                      Class edgeMeshFormat Declaration
\*---------------------------------------------------------------------------*/

class edgeMeshFormat
:
    public edgeMesh
{
protected:

    // Protected Member Functions

        //- Write header information
        static void writeHeader
        (
            Ostream&,
            const pointField&,
            const edgeList&
        );


public:

    // Constructors

        //- Construct from file name
        edgeMeshFormat(const fileName&);

        //- Disallow default bitwise copy construction
        edgeMeshFormat(const edgeMeshFormat&) = delete;


    // Selectors

        //- Read file and return edgeMesh
        static autoPtr<edgeMesh> New(const fileName& name)
        {
            return autoPtr<edgeMesh>
            (
                new edgeMeshFormat(name)
            );
        }


    //- Destructor
    virtual ~edgeMeshFormat()
    {}


    // Member Functions

        //- Read edgeMesh components from stream
        static bool read
        (
            Istream&,
            pointField&,
            edgeList&
        );

        //- Write edgeMesh components to stream
        static Ostream& write
        (
            Ostream&,
            const pointField&,
            const edgeList&,
            bool writeComment =  true
        );

        //- Write edgeMesh with a mimicked IOobject header
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
        void operator=(const edgeMeshFormat&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
