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
    Foam::fvMeshMover

Description
    Abstract base class for fvMesh movers.

    These classes move the mesh points, update the cell volumes and generate
    the corresponding mesh fluxes without any topology change.

SourceFiles
    fvMeshMover.C
    fvMeshMoverNew.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshMover_H
#define fvMeshMover_H

#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fvMeshMover Declaration
\*---------------------------------------------------------------------------*/

class fvMeshMover
{
    // Private Data

        //- Non-const fvMesh reference to allow update
        fvMesh& mesh_;

        //- The dynamicMeshDict
        dictionary dynamicMeshDict_;


public:

    //- Runtime type information
    TypeName("fvMeshMover");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fvMeshMover,
            fvMesh,
            (fvMesh& mesh),
            (mesh)
        );


    // Constructors

        //- Construct from fvMesh
        explicit fvMeshMover(fvMesh&);

        //- Disallow default bitwise copy construction
        fvMeshMover(const fvMeshMover&) = delete;


    // Selectors

        //- Select, construct and return the fvMeshMover.
        //  If the constant/dynamicMeshDict does not exist,
        //  a dummy fvMeshMovers::none object is returned.
        static autoPtr<fvMeshMover> New(fvMesh&);


    //- Destructor
    virtual ~fvMeshMover();


    // Member Functions

        //- Access the fvMesh
        fvMesh& mesh()
        {
            return mesh_;
        }

        //- Return the fvMesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return the dynamicMeshDict/mover dict
        const dictionary& dict() const
        {
            return dynamicMeshDict_.subDict("mover");
        }

        //- Is mesh dynamic, i.e. might it change?
        //  Defaults to true, set to false in the fvMeshMovers::none
        virtual bool dynamic() const
        {
            return true;
        }

        //- Update the mesh for both mesh motion and topology change
        virtual bool update() = 0;

        //-  Update local data for topology changes
        virtual void topoChange(const polyTopoChangeMap&) = 0;

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&) = 0;

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap& mdpm) = 0;

        //- Write the mover state
        virtual bool write(const bool write = true) const
        {
            return true;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvMeshMover&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
