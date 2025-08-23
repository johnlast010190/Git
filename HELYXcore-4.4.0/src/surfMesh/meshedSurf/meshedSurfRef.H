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

Class
    Foam::meshedSurfRef

Description
    Implements a meshed surface by referencing existing faces and points.

\*---------------------------------------------------------------------------*/

#ifndef meshedSurfRef_H
#define meshedSurfRef_H

#include "meshedSurf/meshedSurf.H"
#include "primitives/ints/lists/labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class meshedSurfRef Declaration
\*---------------------------------------------------------------------------*/

class meshedSurfRef
:
    public meshedSurf
{
    const pointField& points_;
    const faceList&   faces_;
    const labelList&  zoneIds_;


    // Private Member Functions

        //- Disallow construct as copy
        meshedSurfRef(const meshedSurfRef&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const meshedSurfRef&) = delete;

public:

    // Public Member Functions

    // Constructors

        //- Construct from components
        meshedSurfRef
        (
            const pointField& pts,
            const faceList& faces,
            const labelList& ids = Foam::emptyLabelList
        )
        :
            points_(pts),
            faces_(faces),
            zoneIds_(ids)
        {}


    //- Destructor
    virtual ~meshedSurfRef()
    {}


    // Access Member Functions

        //- Const access to (global) points used for the surface
        virtual const pointField& points() const
        {
            return points_;
        }

        //- Const access to the surface faces
        virtual const faceList& faces() const
        {
            return faces_;
        }

        //- Const access to per-face zone/region information
        virtual const labelList& zoneIds() const
        {
            return zoneIds_;
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
