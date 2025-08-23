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
    Foam::fvMeshMovers:motionSolver

Description
    The motionSolver fvMeshMover.

SourceFiles
    fvMeshMoversMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshMoversMotionSolver_H
#define fvMeshMoversMotionSolver_H

#include "fvMesh/fvMeshMovers/fvMeshMover/fvMeshMover.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class motionSolver;

namespace fvMeshMovers
{

/*---------------------------------------------------------------------------*\
                  Class fvMeshMovers::motionSolver Declaration
\*---------------------------------------------------------------------------*/

class motionSolver
:
    public fvMeshMover
{
    // Private Data

        autoPtr<Foam::motionSolver> motionPtr_;


public:

    //- Runtime type information
    TypeName("motionSolver");


    // Constructors

        //- Construct from fvMesh
        motionSolver(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        motionSolver(const motionSolver&) = delete;


    //- Destructor
    virtual ~motionSolver();


    // Member Functions

        //- Return the motionSolver
        const Foam::motionSolver& motion() const;

        //- Update the mesh for both mesh motion and topology change
        virtual bool update() override;

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&) override;

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&) override;

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&) override;

        //- Write the motion solver state
        virtual bool write(const bool write = true) const override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const motionSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshMovers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
