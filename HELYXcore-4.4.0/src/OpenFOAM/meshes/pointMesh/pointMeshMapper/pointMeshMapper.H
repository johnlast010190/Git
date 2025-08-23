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
    (c) 2024 Engys Ltd.

Class
    Foam::pointMeshMapper

Description
    Class holds all the necessary information for mapping fields associated
    with pointMesh.

SourceFiles
    pointMeshMapper.C

\*---------------------------------------------------------------------------*/

#ifndef pointMeshMapper_H
#define pointMeshMapper_H

#include "meshes/pointMesh/pointMeshMapper/pointMapper.H"
#include "meshes/pointMesh/pointMeshMapper/pointBoundaryMeshMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class pointMesh;
class polyTopoChangeMap;

/*---------------------------------------------------------------------------*\
                        Class pointMeshMapper Declaration
\*---------------------------------------------------------------------------*/

class pointMeshMapper
{
    // Private data

        //- Reference to mesh
        const pointMesh& mesh_;

        //- Point mapper
        pointMapper pointMap_;

        //- Boundary mapper
        pointBoundaryMeshMapper boundaryMap_;


public:

    // Constructors

        //- Construct from pointMesh
        pointMeshMapper(const pointMesh& mesh, const polyTopoChangeMap& map)
        :
            mesh_(mesh),
            pointMap_(mesh, map),
            boundaryMap_(mesh, pointMap_, map)
        {}

        //- Disallow default bitwise copy construction
        pointMeshMapper(const pointMeshMapper&) = delete;


    // Member Functions


        //- Return reference to mesh fields belong to
        const pointMesh& mesh() const
        {
            return mesh_;
        }

        //- Return reference to objectRegistry storing fields. Can be
        //  removed once fields stored on pointMesh.
        const objectRegistry& thisDb() const
        {
            return mesh_();
        }

        //- Return point mapper
        const morphFieldMapper& pointMap() const
        {
            return pointMap_;
        }

        //- Return boundary mapper
        const pointBoundaryMeshMapper& boundaryMap() const
        {
            return boundaryMap_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pointMeshMapper&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
