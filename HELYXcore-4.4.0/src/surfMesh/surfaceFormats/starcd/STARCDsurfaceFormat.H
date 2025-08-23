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
    Foam::fileFormats::STARCDsurfaceFormat

Description
    Read/write the surface shells from pro-STAR vrt/cel files.

Note
    Uses the extension \a .inp (input) to denote the format.

See also
    Foam::meshReaders::STARCD

SourceFiles
    STARCDsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef STARCDsurfaceFormat_H
#define STARCDsurfaceFormat_H

#include "MeshedSurface/MeshedSurface.H"
#include "MeshedSurfaceProxy/MeshedSurfaceProxy.H"
#include "UnsortedMeshedSurface/UnsortedMeshedSurface.H"
#include "surfaceFormats/starcd/STARCDsurfaceFormatCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                     Class STARCDsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class STARCDsurfaceFormat
:
    public MeshedSurface<Face>,
    public STARCDsurfaceFormatCore
{
    // Private Member Functions

        static inline void writeShell
        (
            Ostream&,
            const Face&,
            const label cellId,
            const label cellTableId
        );


public:

    // Constructors

        //- Construct from file name
        STARCDsurfaceFormat(const fileName&);

        //- Disallow default bitwise copy construction
        STARCDsurfaceFormat(const STARCDsurfaceFormat<Face>&) = delete;


    // Selectors

        //- Read file and return surface
        static autoPtr<MeshedSurface<Face>> New(const fileName& name)
        {
            return autoPtr<MeshedSurface<Face>>
            (
                new STARCDsurfaceFormat<Face>(name)
            );
        }


    //- Destructor
    virtual ~STARCDsurfaceFormat()
    {}


    // Member Functions

        //- Write surface mesh components by proxy
        static void write
        (
            const fileName&,
            const MeshedSurfaceProxy<Face>&
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
        void operator=(const STARCDsurfaceFormat<Face>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "surfaceFormats/starcd/STARCDsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
