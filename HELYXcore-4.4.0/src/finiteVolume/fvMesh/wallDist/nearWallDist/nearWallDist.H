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
    (c) 2023-2024 Engys Ltd.

Class
    Foam::nearWallDist

Description
    Distance calculation for cells with face on a wall.
    Searches pointNeighbours to find closest.

SourceFiles
    nearWallDist.C

\*---------------------------------------------------------------------------*/

#ifndef nearWallDist_H
#define nearWallDist_H

#include "meshes/MeshObject/MeshObject.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                        Class nearWallDist Declaration
\*---------------------------------------------------------------------------*/

class nearWallDist
:
    public MeshObject<fvMesh, UpdateableMeshObject, nearWallDist>
{
    // Private data

        //- Near-well distance-to-wall boundary field
        volScalarField::Boundary y_;


    // Private Member Functions

        //- Re-size y_ following a change to the boundary mesh
        void resize();

        //- Re-calculate y_
        void correct();


public:

    // Declare name of the class and its debug switch
    ClassName("nearWallDist");


    // Constructors

        //- Construct from components
        nearWallDist(const fvMesh& mesh);

        //- Disallow default bitwise copy construction
        nearWallDist(const nearWallDist&) = delete;


    //- Destructor
    virtual ~nearWallDist();


    // Member Functions

        //- Access to the near-wall distances
        const volScalarField::Boundary& y() const
        {
            return y_;
        }

        //- Update the y-field when the mesh moves
        virtual bool movePoints();

        //- Update the y-field when the mesh changes
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const nearWallDist&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
