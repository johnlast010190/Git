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
    (c) 2025 Engys Ltd.
    (c) 2011-2016 OpenFOAM Foundation

Class
    Foam::deformingBodyFvMotionSolver

Group
    grpMeshMotionSolvers

Description
    Deforming body mesh motion solver for an fvMesh.

    Applies motion description extrapolation of movement as function of
    distance to the object surface. The algorithm is an adaptation of
    rigidBodyMeshMotion.

SourceFiles
    deformingBodyFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef deformingBodyFvMotionSolver_H
#define deformingBodyFvMotionSolver_H

#include "motionSolvers/displacement/displacement/displacementMotionSolver.H"
#include "fvMotionSolvers/fvMotionSolver/fvMotionSolver.H"
#include "meshes/polyMesh/polyTopoChangeMap/polyTopoChangeMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class deformingBodyFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class deformingBodyFvMotionSolver
:
    public displacementMotionSolver,
    public fvMotionSolver
{
    //- Class containing the patches and point motion weighting for each body
    class bodyMesh
    {
        //- Name of the body
        const word name_;

        //- List of mesh patches associated with this body
        const wordReList patches_;

        //- Patches to integrate forces
        const labelHashSet patchSet_;

        //- Inner morphing distance (limit of solid-body region)
        const scalar di_;

        //- Outer morphing distance (limit of linear interpolation region)
        const scalar do_;

        //- Current interpolation weight
        //  (1 at patches and within di_, 0 at do_ and beyond)
        pointScalarField weight_;


    public:

        friend class deformingBodyFvMotionSolver;

        bodyMesh
        (
            const polyMesh& mesh,
            const word& name,
            const dictionary& dict
        );

        template<class Type>
        inline Type weight(const Type& pDist) const;
    };


protected:
    // Protected data

        //- List of the bodyMeshes containing the patches and point motion
        //  weighting for each body
        PtrList<bodyMesh> bodyMeshes_;

        //- Optional cell-centre motion field
        autoPtr<volVectorField> cellDisplacement_;

        //- Optional write cellDisplacement field and weights
        bool writeFields_;

        //- Return the list of weights for each body to pointi
        //   w is used as the work-space for the weights and returned
        List<scalar>& weights(const label pointi, List<scalar>& w) const;

        //- Return the analytical displacement for body N
        tmp<pointField> analyticalFieldDisplacement(const bodyMesh&);

        //- Initialise the bodyMeshes
        void initBodies(const polyMesh&, const dictionary&);

        //- Initialise optional fields to write
        void initFields();


public:

    //- Runtime type information
    TypeName("deformingBodyMotion");


    // Constructors

        //- Construct from polyMesh and dictionary
        deformingBodyFvMotionSolver
        (
            const polyMesh&,
            const dictionary&
        );

        //- Construct from components
        deformingBodyFvMotionSolver
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const pointVectorField& pointDisplacement,
            const pointVectorField& points0
        );

        //- Disallow default bitwise copy construction
        deformingBodyFvMotionSolver
        (
            const deformingBodyFvMotionSolver&
        ) = delete;


    //- Destructor
    ~deformingBodyFvMotionSolver();


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const deformingBodyFvMotionSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
