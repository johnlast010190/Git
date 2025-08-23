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
    (c) 2021-2022 OpenFOAM Foundation
    (c) 2023 Engys Ltd.

Class
    Foam::fvMeshMovers::none

Description
    Dummy fvMeshMover which does not move the mesh points.

SourceFiles
    none.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshMoversNone_H
#define fvMeshMoversNone_H

#include "fvMesh/fvMeshMovers/fvMeshMover/fvMeshMover.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshMovers
{

/*---------------------------------------------------------------------------*\
                      Class fvMeshMovers::none Declaration
\*---------------------------------------------------------------------------*/

class none
:
    public fvMeshMover
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from fvMesh
        none(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        none(const none&) = delete;


    //- Destructor
    virtual ~none();


    // Member Functions

        //- The mesh is static and does not change, so return false
        virtual bool dynamic() const override
        {
            return false;
        }

        //- Dummy update function which does not change the mesh
        virtual bool update() override
        {
            return false;
        }

        //-  Update local data for topology changes
        virtual void topoChange(const polyTopoChangeMap&) override
        {}

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&) override
        {}

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&) override
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const none&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshMovers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
