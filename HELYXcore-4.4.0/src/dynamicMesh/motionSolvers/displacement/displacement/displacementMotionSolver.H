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
    (c) 2024 Engys Ltd.

Class
    Foam::displacementMotionSolver

Description
    Virtual base class for displacement motion solver

    The boundary displacement is set as a boundary condition
    on the pointDisplacement pointVectorField.

SourceFiles
    displacementMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef displacementMotionSolver_H
#define displacementMotionSolver_H

#include "motionSolvers/displacement/points/pointsMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class displacementMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class displacementMotionSolver
:
    public pointsMotionSolver
{
protected:

    // Protected data

        //- Point motion field
        mutable pointVectorField pointDisplacement_;


public:

    //- Runtime type information
    TypeName("displacementMotionSolver");


    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        displacementMotionSolver,
        displacement,
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const pointVectorField& pointDisplacement,
            const pointVectorField& points0
        ),
        (mesh, dict, pointDisplacement, points0)
    );


    // Selectors

        //- Select constructed from polyMesh, dictionary and components.
        //  If dictionary was registered this will 'steal' that registration.
        static autoPtr<displacementMotionSolver> New
        (
            const word& solverTypeName,
            const polyMesh&,
            const dictionary&,
            const pointVectorField& pointDisplacement,
            const pointVectorField& points0
        );


    // Constructors

        //- Construct from mesh and dictionary
        displacementMotionSolver
        (
            const polyMesh&,
            const dictionary&,
            const word& type
        );

        //- Construct from components
        displacementMotionSolver
        (
            const polyMesh&,
            const dictionary&,
            const pointVectorField& pointDisplacement,
            const pointVectorField& points0,
            const word& type
        );

        //- Disallow default bitwise copy construction
        displacementMotionSolver(const displacementMotionSolver&) = delete;


    //- Destructor
    virtual ~displacementMotionSolver();


    // Member Functions

        //- Return reference to the point motion displacement field
        pointVectorField& pointDisplacement()
        {
            return pointDisplacement_;
        }

        //- Return const reference to the point motion displacement field
        const pointVectorField& pointDisplacement() const
        {
            return pointDisplacement_;
        }

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const displacementMotionSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
