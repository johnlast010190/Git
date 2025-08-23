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
    Foam::fileFormats::OBJsurfaceFormat

Description
    Provide a means of reading/writing Alias/Wavefront OBJ format.

    Does not handle negative face indices.

SourceFiles
    OBJsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef OBJsurfaceFormat_H
#define OBJsurfaceFormat_H

#include "MeshedSurface/MeshedSurface.H"
#include "MeshedSurfaceProxy/MeshedSurfaceProxy.H"
#include "UnsortedMeshedSurface/UnsortedMeshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                      Class OBJsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class OBJsurfaceFormat
:
    public MeshedSurface<Face>
{
public:

    // Constructors

        //- Construct from file name
        OBJsurfaceFormat(const fileName&);

        //- Disallow default bitwise copy construction
        OBJsurfaceFormat(const OBJsurfaceFormat<Face>&) = delete;


    // Selectors

        //- Read file and return surface
        static autoPtr<MeshedSurface<Face>> New(const fileName& name)
        {
            return autoPtr<MeshedSurface<Face>>
            (
                new OBJsurfaceFormat<Face>(name)
            );
        }


    //- Destructor
    virtual ~OBJsurfaceFormat()
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

        //- Write object file
        virtual void write(const fileName& name) const
        {
            write(name, MeshedSurfaceProxy<Face>(*this));
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const OBJsurfaceFormat<Face>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "surfaceFormats/obj/OBJsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
