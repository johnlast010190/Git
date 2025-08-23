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
    Foam::fileFormats::VTKedgeFormat

Description
    Provide a means of writing VTK legacy format.

SourceFiles
    VTKedgeFormat.C

\*---------------------------------------------------------------------------*/

#ifndef VTKedgeFormat_H
#define VTKedgeFormat_H

#include "edgeMesh/edgeMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                     Class VTKedgeFormat Declaration
\*---------------------------------------------------------------------------*/

class VTKedgeFormat
:
    public edgeMesh
{
protected:

    // Protected Member Functions

        //- Write header information with points
        static void writeHeader
        (
            Ostream&,
            const pointField&
        );

        //- Write edges
        static void writeEdges(Ostream&, const UList<edge>&);


public:


    // Constructors

        //- Construct from file name
        VTKedgeFormat(const fileName&);

        //- Disallow default bitwise copy construction
        VTKedgeFormat(const VTKedgeFormat&) = delete;


    // Selectors

        //- Read file and return surface
        static autoPtr<edgeMesh> New(const fileName& name)
        {
            return autoPtr<edgeMesh>
            (
                new VTKedgeFormat(name)
            );
        }


    //- Destructor
    virtual ~VTKedgeFormat()
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
        void operator=(const VTKedgeFormat&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
