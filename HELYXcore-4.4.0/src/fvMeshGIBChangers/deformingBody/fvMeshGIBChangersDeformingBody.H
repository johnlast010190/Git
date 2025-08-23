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
    (c) 2024 Engys Ltd.

Class
    Foam::fvMeshGIBChangers::deformingBody

Description
    Dynamic GIB mesh changer for deforming body motion.

SourceFiles
    fvMeshGIBChangersDeformingBody.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshGIBChangersDeformingBody_H
#define fvMeshGIBChangersDeformingBody_H

#include "base/fvMeshGIBChangersBase.H"
#include "deformingBody/motionFunctions/motionFunction/deformingBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshGIBChangers
{

/*---------------------------------------------------------------------------*\
               Class fvMeshGIBChangers::deformingBody Declaration
\*---------------------------------------------------------------------------*/

class deformingBody
:
    public fvMeshGIBChangersBase
{
    // Private Member Functions

        virtual void initialization();

        virtual void makeRecAllPoints0() const override;


protected:

    // Protected Data

        //- stl interface stl location
        fileName postProFolder_;

        //- Interface output switch
        Switch writeInterface_;

        //- Interface properties
        autoPtr<motionFunction> DBMFPtr_;


    // Protected Member Functions

        virtual bool includeWalls() const override
        {
            return false;
        }

        //- Do tangential boundary motion
        virtual bool tangentialBoundaryMotion() const override
        {
            return true;
        }

        //- Move old GIB
        virtual tmp<vectorField> movePolyPatch
        (
            primitivePatch&,
            const vectorField&
        ) override;

        void writeInterface(const indirectPolyPatch&);

        //- postPro folder
        void postPro();


public:

    //- Runtime type information
    TypeName("deformingBody");


    // Constructors

        //- Construct from fvMesh
        explicit deformingBody(fvMesh&);

        //- Construct from fvMesh and type
        deformingBody(fvMesh&, const word&);

        //- Construct from fvMesh and dictionary
        deformingBody(fvMesh&, const dictionary);

        //- Disallow default bitwise copy construction
        deformingBody(const deformingBody&) = delete;


    //- Destructor
    ~deformingBody();


    // Member Functions

        virtual void updateInit(const word&) override;

        //- Update the mesh for both mesh motion and topology change
        virtual bool update() override;

        //- Do tangential boundary motion
        virtual bool tangentialBoundaryMotion()
        {
            return true;
        }

        //- Correct velocity for GIB
        virtual tmp<vectorField> velocityCorrect
        (
            const vectorField&
        ) const override;

        //- Calculate the old boundary location
        virtual tmp<vectorField> oldBoundaryLocation() const override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const deformingBody&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshGIBChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
