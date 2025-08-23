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
    (c) 2024-2025 Engys Ltd.

Class
    Foam::fvMeshGIBChangers::adjoint

Description
    Dynamic GIB mesh changer for adjoint method.

SourceFiles
    fvMeshGIBChangersAdjoint.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshGIBChangersAdjoint_H
#define fvMeshGIBChangersAdjoint_H

#include "deformingBody/fvMeshGIBChangersDeformingBody.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshGIBChangers
{

/*---------------------------------------------------------------------------*\
                  Class fvMeshGIBChangers::adjoint Declaration
\*---------------------------------------------------------------------------*/

class adjoint
:
    public deformingBody
{
    // Private Data

        //- Reference to adjoint properties
        const IOdictionary* adjPropertiesPtr_;


    // Private Member Functions

        //- Courant-based scaling based on a max and a mean courant
        void courantScaling(scalarField&);

        //- Point curvature calculation
        virtual scalarField calcCurvature();

        //- If new points of the poly-patch are moving out of the base mesh
        //  bounds, find the intersections at the base mesh boundary and
        //  then snap there (!!! round errors careful).
        virtual void constraintPolyPointsOutMotion(vectorField&);

        virtual void nearBoundaryIntersectionsChecking(vectorField&);

        virtual void checkConcaveBoundaryPatchPoints(vectorField&);

        //- Positive point sense at boundary points, which are on the GIB
        //  with positive point sensitivity fixed to zero.
        virtual void filterBoundaryPoints
        (
            scalarField&,
            const indirectPolyPatch&
        );

        //- Mark the boundary faces of the GIB that have not moved and
        //  the sensitivity want to move outwards.
        //  These faces are constrained inside the smoother.
        //  COUPLED needs correction.
        virtual boolList findConstraintFaces();

        //- Fix the points of the constraint patches (if moved).
        //  Only inlet outlet currently supported.
        virtual void fixConstraintPatches(pointField&);

        //- If processor patch point changed location, then sync the other side.
        virtual void syncProcBoundaryPoints(pointField&, const pointField&);


protected:

        virtual bool includeWalls() const override
        {
            return true;
        }

        //- Do tangential boundary motion
        virtual bool tangentialBoundaryMotion() const override
        {
            return false;
        }

        //- Move old GIB based on the curvature and the sensitivities
        virtual tmp<vectorField> computeNewPoints
        (
            primitivePatch&,
            const scalarField&
        );


public:

    //- Runtime type information
    TypeName("adjoint");


    // Constructors

        //- Construct from fvMesh
        explicit adjoint(fvMesh&);

        //- Construct from fvMesh and dictionary
        adjoint(fvMesh&, const dictionary);

        //- Disallow default bitwise copy construction
        adjoint(const adjoint&) = delete;


    //- Destructor
    ~adjoint();


    // Member Functions

        //- faceZone initialization using the stl information
        virtual void updateInit();

        //- Update the gib
        virtual bool update() override;

        //- Reference to the adjoint properties
        const IOdictionary& adjProperties();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const adjoint&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshGIBChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
