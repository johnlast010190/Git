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
    Foam::meshedSurf

Description
    Abstract definition of a meshed surface defined by faces and points.

\*---------------------------------------------------------------------------*/

#ifndef meshedSurf_H
#define meshedSurf_H

#include "meshes/primitiveShapes/point/pointField.H"
#include "meshes/meshShapes/face/faceList.H"
#include "containers/Lists/ListOps/ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class meshedSurf Declaration
\*---------------------------------------------------------------------------*/

class meshedSurf
{
public:

    // Constructors

        //- Construct null
        meshedSurf()
        {}


    //- Destructor
    virtual ~meshedSurf()
    {}


    // Access Member Functions

        //- Const access to (global) points used for the surface
        virtual const pointField& points() const = 0;

        //- Const access to the surface faces
        virtual const faceList& faces() const = 0;

        //- Const access to per-face zone/region information
        virtual const labelList& zoneIds() const = 0;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
