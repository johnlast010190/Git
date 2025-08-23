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
    (c) 2010-2012 Engys Ltd.
    (c) 1991-2007 OpenCFD Ltd.

Class
    ensightWriter

Description

SourceFiles
    ensightWriter.C

\*---------------------------------------------------------------------------*/

#ifndef ensightWriter_H
#define ensightWriter_H

#include "triSurfaceWriter/triSurfaceWriter.H"
#include "primitives/Scalar/floatScalar/floatScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class triSurface;

/*---------------------------------------------------------------------------*\
                           Class ensightWriter Declaration
\*---------------------------------------------------------------------------*/

class ensightWriter
:
    public triSurfaceWriter
{

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("ensight");


    // Constructors

        //- Construct null
        ensightWriter();

    // Destructor

        virtual ~ensightWriter();


    // Member Functions

        //- Write surface + data to file.
        virtual void write
        (
            const fileName& path,           // <root>/<case>/sampleSurfaces
            const fileName& timeDir,        // time name
            const fileName& surfaceName,    // name of surface
            const triSurface& surface,
            PtrList<triSurfaceScalarField>& scalarFields,
            PtrList<triSurfaceVectorField>& vectorFields,
            PtrList<triSurfaceTensorField>& tensorFields
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
