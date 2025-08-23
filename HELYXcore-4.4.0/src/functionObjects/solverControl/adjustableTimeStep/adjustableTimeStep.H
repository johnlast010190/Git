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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2010-2024 Engys Ltd

Class
    adjustableTimeStep

Group
    grpSolverControlFunctionObjects

Description
    Run time control of the time step

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    adjustableTimeStep.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_adjustableTimeStep_H
#define functionObjects_adjustableTimeStep_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class adjustableTimeStep Declaration
\*---------------------------------------------------------------------------*/

class adjustableTimeStep
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        //- input
        scalar targetCo_;
        scalar maxDeltaT_;
        scalar minDeltaT_;

        word CoType_;
        word phiName_;

        //- Result values
        scalar deltaT_;

    // Private Member Functions

        scalar getCoNum();


public:

    //- Runtime type information
    TypeName("adjustableTimeStep");


    // Constructors

        //- Construct for given objectRegistry and dictionary
        adjustableTimeStep
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        adjustableTimeStep(const adjustableTimeStep&) = delete;


    //- Destructor
    virtual ~adjustableTimeStep();


    // Member Functions

        //- Update solver
        void calculate();

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Read the parameters
        virtual bool read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const adjustableTimeStep&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
