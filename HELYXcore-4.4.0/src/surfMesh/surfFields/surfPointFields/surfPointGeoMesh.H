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
    (c) 2011 OpenFOAM Foundation

Class
    Foam::surfPointGeoMesh

Description
    The surfMesh GeoMesh (for holding fields).

    Similar to surfGeoMesh, but refers to the surface points.

\*---------------------------------------------------------------------------*/

#ifndef surfPointGeoMesh_H
#define surfPointGeoMesh_H

#include "meshes/GeoMesh/GeoMesh.H"
#include "surfMesh/surfMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class surfPointGeoMesh Declaration
\*---------------------------------------------------------------------------*/

class surfPointGeoMesh
:
    public GeoMesh<surfMesh>
{

public:

    // Constructors

        //- Construct from surfMesh reference
        explicit surfPointGeoMesh(const surfMesh& mesh)
        :
            GeoMesh<surfMesh>(mesh)
        {}


    // Member Functions

        //- Return size
        static label size(const surfMesh& mesh)
        {
            return mesh.nPoints();
        }

        //- Return size
        label size() const
        {
            return size(mesh_);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
