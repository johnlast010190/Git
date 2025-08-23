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
    (c) 2012-2016 OpenFOAM Foundation
    (c) 2023 FOSS GP
    (c) 2024 Engys Ltd.

Class
    Foam::modesMotionSolver

Description
    Reads in morphing modes and deforms the mesh by multipling them with
    given amplifications

SourceFiles
    modesMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef modesMotionSolver_H
#define modesMotionSolver_H

#include "motionSolvers/motionSolver/motionSolver.H"
#include "fields/GeometricFields/pointFields/pointFields.H"
#include "primitives/strings/lists/wordReList.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class modesMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class modesMotionSolver
:
    public motionSolver
{
protected:

    // Protected data

        //- List of names containing the modes
        wordReList fileNames_;

        //- Amplification of each mode
        PtrList<Function1<scalar>> amplifications_;

        //- The modes
        PtrList<pointVectorField> modes_;

        //- The initial points
        vectorField points0_;

        //- Append to movement to the current points
        bool appendToCurrentPoints_;


private:

    // Private Member Functions

        //- No copy construct
        modesMotionSolver
        (
            const modesMotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const modesMotionSolver&) = delete;

public:

    //- Runtime type information
    TypeName("modesMotionSolver");


    // Constructors

        //- Construct from mesh and dictionary
        modesMotionSolver
        (
            const polyMesh&,
            const dictionary&
        );


    //- Destructor
    virtual ~modesMotionSolver() = default;


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        //  Does nothing
        virtual void solve(){};

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&);

        //-  Update local data for topology changes
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&){}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
