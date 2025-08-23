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
    (c) 2016 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.
    (c) 2023-2024 Engys Ltd.

Class
    Foam::pointsMotionSolver

Description
    Virtual base class for displacement motion solvers, where the point
    motion is relative to a set of fixed points (points0).

SourceFiles
    pointsMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef pointsMotionSolver_H
#define pointsMotionSolver_H

#include "motionSolvers/motionSolver/motionSolver.H"
#include "fields/GeometricFields/pointFields/pointFields.H"
//#include "meshes/primitiveShapes/point/pointIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyTopoChangeMap;

/*---------------------------------------------------------------------------*\
                     Class pointsMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class pointsMotionSolver
:
    public motionSolver
{
protected:

    // Protected data

        //- Starting points
        mutable pointVectorField points0_;


public:

    //- Runtime type information
    TypeName("pointsMotionSolver");


    // Constructors

        //- Construct from mesh and dictionary
        pointsMotionSolver
        (
            const polyMesh&,
            const dictionary&,
            const word& type
        );

        //- Construct from mesh, dictionary and fixed points field
        pointsMotionSolver
        (
            const polyMesh&,
            const dictionary&,
            const pointVectorField& points0,
            const word& type
        );

        //- Disallow default bitwise copy construction
        pointsMotionSolver(const pointsMotionSolver&) = delete;


    //- Destructor
    virtual ~pointsMotionSolver();


    // Member Functions

        static pointVectorField readPoints0(const polyMesh& mesh);

        //- Return reference to the reference field
        pointField& points0()
        {
            return points0_;
        }

        //- Return reference to the reference field
        const pointField& points0() const
        {
            return points0_;
        }

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&);

        //-  Update local data for topology changes
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map.
        //  Resets points0 to the points of the new mesh
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);

        //- Write points0 if the mesh topology changed
        virtual bool write() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pointsMotionSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
