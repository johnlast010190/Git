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
    (c) 2020-2024 Engys Ltd.

Class
    Foam::multiSolidBodyGIBMotionSolver

Description
    Solid-body motion of the mesh specified by a run-time selectable motion
    function, returning STL deformation for use in fvMeshGIBChangersSolidBody
    motion type.

SourceFiles
    multiSolidBodyGIBMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef multiSolidBodyGIBMotionSolver_H
#define multiSolidBodyGIBMotionSolver_H

#include "motionSolvers/motionSolver/motionSolver.H"
#include "motionSolvers/displacement/solidBody/multiSolidBodyMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class multiSolidBodyGIBMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class multiSolidBodyGIBMotionSolver
:
    public motionSolver
{
    // Private data

        //- Specified cellZones
        labelList zoneIDs_;

        //- Points to move per cellZone
        labelListList pointIDs_;

        //- name of the stl geometry
        const word triName_;

        //- Original stl points
        pointField triPoints0_;


public:

    //- Runtime type information
    TypeName("multiSolidBodyGIBMotionSolver");


    // Constructors

        //- Construct from mesh and dictionary
        multiSolidBodyGIBMotionSolver
        (
            const polyMesh&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        multiSolidBodyGIBMotionSolver
        (
            const multiSolidBodyGIBMotionSolver&
        ) = delete;


    //- Destructor
    ~multiSolidBodyGIBMotionSolver();


    // Member Functions

        //- Return point location of the STL obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve()
        {}

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&)
        {}

        //-  Update local data for topology changes
        virtual void topoChange(const polyTopoChangeMap&)
        {}

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&)
        {}

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&)
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multiSolidBodyGIBMotionSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
