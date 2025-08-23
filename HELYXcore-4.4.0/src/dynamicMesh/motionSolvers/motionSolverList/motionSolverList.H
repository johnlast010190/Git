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
    (c) 2019 OpenFOAM Foundation
    (c) 2023-2024 Engys Ltd.

Class
    Foam::motionSolverList

Description
    Motion of the mesh specified as a list of motion solvers.

    The motion solvers are executed in order and the resulting displacements
    accumulated into an overall displacement and the displaced point positions
    returned.

SourceFiles
    motionSolverList.C

\*---------------------------------------------------------------------------*/

#ifndef motionSolverList_H
#define motionSolverList_H

#include "motionSolvers/motionSolver/motionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class motionSolverList Declaration
\*---------------------------------------------------------------------------*/

class motionSolverList
:
    public motionSolver
{
    // Private data

        PtrList<motionSolver> motionSolvers_;


public:

    //- Runtime type information
    TypeName("motionSolverList");


    // Constructors

        //- Construct from mesh and dictionary
        motionSolverList(const polyMesh&, const dictionary&);

        //- Disallow default bitwise copy construction
        motionSolverList(const motionSolverList&) = delete;


    //- Destructor
    virtual ~motionSolverList();


    // Member Functions

        //- Provide current points for motion. Use current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&);

        //-  Update local data for topology changes
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);

        //- Write motion state information for restart
        virtual bool write() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const motionSolverList&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
