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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::fileFormats::STLsurfaceFormat

Description
    Provide a means of reading/writing STL files (ASCII and BINARY).

Note
    For efficiency, the zones are sorted before creating the faces.
    The class is thus derived from MeshedSurface.

SourceFiles
    STLsurfaceFormat.C
    STLsurfaceFormatASCII.L

\*---------------------------------------------------------------------------*/

#ifndef STLsurfaceFormat_H
#define STLsurfaceFormat_H

#include "stl/STLReader.H"
#include "MeshedSurface/MeshedSurface.H"
#include "MeshedSurfaceProxy/MeshedSurfaceProxy.H"
#include "UnsortedMeshedSurface/UnsortedMeshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                     Class STLsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class STLsurfaceFormat
:
    public MeshedSurface<Face>,
    public STLCore
{
    // Private Member Functions

        //- Write Face (ASCII)
        static inline void writeShell
        (
            Ostream&,
            const pointField&,
            const Face&
        );

        //- Write Face (BINARY)
        static inline void writeShell
        (
            ostream&,
            const pointField&,
            const Face&,
            const label zoneI
        );


public:

    // Constructors

        //- Construct from file name
        STLsurfaceFormat(const fileName&);

        //- Disallow default bitwise copy construction
        STLsurfaceFormat(const STLsurfaceFormat<Face>&) = delete;


    // Selectors

        //- Read file and return surface
        static autoPtr<MeshedSurface<Face>> New(const fileName& name)
        {
            return autoPtr<MeshedSurface<Face>>
            (
                new STLsurfaceFormat<Face>(name)
            );
        }


    //- Destructor
    virtual ~STLsurfaceFormat()
    {}


    // Member Functions

        //- Write surface mesh components by proxy (as ASCII)
        static void writeAscii
        (
            const fileName&,
            const MeshedSurfaceProxy<Face>&
        );

        //- Write surface mesh components by proxy (as BINARY)
        static void writeBinary
        (
            const fileName&,
            const MeshedSurfaceProxy<Face>&
        );

        //- Write surface mesh components by proxy
        //  as ASCII or BINARY, depending on the extension
        static void write
        (
            const fileName&,
            const MeshedSurfaceProxy<Face>&
        );

        //- Write surface mesh components by proxy
        //  as ASCII or BINARY or dependent on the extension
        static void write
        (
            const fileName&,
            const MeshedSurfaceProxy<Face>&,
            const STLFormat&
        );

        //- Write UnsortedMeshedSurface (as ASCII) sorted by zone
        static void writeAscii
        (
            const fileName&,
            const UnsortedMeshedSurface<Face>&
        );

        //- Write UnsortedMeshedSurface (as BINARY) unsorted by zone
        static void writeBinary
        (
            const fileName&,
            const UnsortedMeshedSurface<Face>&
        );

        //- Write UnsortedMeshedSurface
        //  as ASCII or BINARY, depending on the extension
        static void write
        (
            const fileName&,
            const UnsortedMeshedSurface<Face>&
        );

        //- Write UnsortedMeshedSurface
        //  as ASCII or BINARY or dependent on the extension
        static void write
        (
            const fileName&,
            const UnsortedMeshedSurface<Face>&,
            const STLFormat&
        );

        //- Read from file
        virtual bool read(const fileName&);

        //- Write object
        virtual void write(const fileName& name) const
        {
            write(name, MeshedSurfaceProxy<Face>(*this));
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const STLsurfaceFormat<Face>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "surfaceFormats/stl/STLsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
