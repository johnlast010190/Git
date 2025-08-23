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
    (c) 2023-2024 Engys Ltd.

Class
    Foam::dynamicFvMesh

Description
    Abstract base class for geometry and/or topology changing fvMesh.

SourceFiles
    dynamicFvMesh.C
    dynamicFvMeshNew.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicFvMesh_H
#define dynamicFvMesh_H

#include "fvMesh/fvMesh.H"
#include "memory/autoPtr/autoPtr.H"
#include "nonConformal/polyPatches/nonConformal/nonConformalPolyPatch.H"

#define MOVING_MESH

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dynamicFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicFvMesh
:
    public fvMesh
{
    // Private data

        //- The dynamicMeshDict
        dictionary dynamicMeshDict_;


    // Private Member Functions

        //- Return the dynamicMeshDict IOobject
        static IOobject dynamicMeshDictIOobject(const IOobject& io);


public:

    //- Runtime type information
    TypeName("dynamicFvMesh");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            dynamicFvMesh,
            IOobject,
            (const IOobject& io),
            (io)
        );


    // Constructors

        //- Construct from objectRegistry, and read/write options
        explicit dynamicFvMesh(const IOobject& io);

        //- Construct from components without boundary.
        //  Boundary is added using addFvPatches() member function
        dynamicFvMesh
        (
            const IOobject& io,
            pointField&& points,
            faceList&& faces,
            labelList&& allOwner,
            labelList&& allNeighbour,
            const bool syncPar = true
        );

        //- Construct without boundary from cells rather than owner/neighbour.
        //  Boundary is added using addPatches() member function
        dynamicFvMesh
        (
            const IOobject& io,
            pointField&& points,
            faceList&& faces,
            cellList&& cells,
            const bool syncPar = true
        );

        //- Disallow default bitwise copy construction
        dynamicFvMesh(const dynamicFvMesh&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<dynamicFvMesh> New(const IOobject& io);


    //- Destructor
    virtual ~dynamicFvMesh();


    // Member Functions

        //- Return the dynamicMeshDict
        const dictionary& dynamicMeshDict() const
        {
            return dynamicMeshDict_;
        }

        //- Update the mesh for both topology change and mesh motion
        virtual bool update() = 0;

        //- Update the mesh to the provided configuration
        virtual bool update(const word&);

        //- needed for the initialization of the GIB classes
        virtual void updateInit(){};

        //- Update the mesh for topology change only
        virtual bool topoChange()
        {
            return false;
        }

        //- Update the mesh for topology change only
        virtual void topoChange(const polyTopoChangeMap& map)
        {
            fvMesh::topoChange(map);
        }

        //- Move the mesh
        virtual bool move()
        {
            return update();
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const dynamicFvMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
