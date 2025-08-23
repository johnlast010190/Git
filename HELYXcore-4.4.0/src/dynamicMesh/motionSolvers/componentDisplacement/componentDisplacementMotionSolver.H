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
    Foam::componentDisplacementMotionSolver

Description
    Virtual base class for displacement motion solver

    The boundary displacement is set as a boundary condition
    on the pointDisplacementX pointScalarField.

SourceFiles
    componentDisplacementMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef componentDisplacementMotionSolver_H
#define componentDisplacementMotionSolver_H

#include "motionSolvers/motionSolver/motionSolver.H"
#include "fields/GeometricFields/pointFields/pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyTopoChangeMap;

/*---------------------------------------------------------------------------*\
                   Class componentDisplacementMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class componentDisplacementMotionSolver
:
    public motionSolver
{
protected:

    // Protected data

        //- The component name to solve for
        word cmptName_;

        //- The component to solve for
        direction cmpt_;

        //- Reference point field for this component
        scalarField points0_;

        //- Point motion field
        mutable pointScalarField pointDisplacement_;


private:

    // Private Member Functions

        //- Return the component corresponding to the given component name
        direction cmpt(const word& cmptName) const;


public:

    //- Runtime type information
    TypeName("componentDisplacementMotionSolver");


    // Constructors

        //- Construct from polyMesh and dictionary and type
        componentDisplacementMotionSolver
        (
            const polyMesh&,
            const dictionary&,
            const word& type
        );

        //- Disallow default bitwise copy construction
        componentDisplacementMotionSolver
        (
            const componentDisplacementMotionSolver&
        ) = delete;


    //- Destructor
    virtual ~componentDisplacementMotionSolver();


    // Member Functions

        //- Return reference to the reference field
        scalarField& points0()
        {
            return points0_;
        }

        //- Return reference to the reference field
        const scalarField& points0() const
        {
            return points0_;
        }

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&);

        //-  Update local data for topology changes
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const componentDisplacementMotionSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
