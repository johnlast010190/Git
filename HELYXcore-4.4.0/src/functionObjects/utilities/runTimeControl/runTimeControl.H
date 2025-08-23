/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
|    o     o     |  ENGYS Ltd. <http://engys.com/>
|       o        |
\*---------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    (c) 2016 OpenCFD Ltd.
    (c) 2015-2016 OpenFOAM Foundation

Class
    Foam::functionObjects::runTimeControl

Group
    grpUtilitiesFunctionObjects

Description
    Controls when the calculation is terminated based on satisfying
    user-specified conditions.

    Optionally specify a number of write steps before the calculation is
    terminated.  Here, a write is performed each time that all conditions are
    satisfied.

SourceFiles
    runTimeControl.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimeControl_H
#define functionObjects_runTimeControl_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "containers/HashTables/Map/Map.H"
#include "primitives/enums/Enum.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimeControls
{

// Forward declaration of classes
class runTimeCondition;

/*---------------------------------------------------------------------------*\
                       Class runTimeControl Declaration
\*---------------------------------------------------------------------------*/

class runTimeControl
:
    public fvMeshFunctionObject
{
public:

    // Public enumerations

        enum class satisfiedAction
        {
            ABORT,         //!< "abort" - write and emit a FatalError
            END,           //!< "end" - write and exit cleanly
            SET_TRIGGER,   //!< "setTrigger" - trigger another condition
            MODIFY_SOLVER  //!< "modifySolver" - modify solver settings
        };

        static Enum<satisfiedAction> satisfiedActionNames;

        enum logicalOps
        {
            logicalAND,
            logicalOR
        };

private:

    // Private data

        //- List of conditions to satisfy
        PtrList<runTimeCondition> conditions_;

        //- Map to define group IDs
        Map<label> groupMap_;
        labelList groupList_;

        //- Number of write steps before exiting
        label nWriteStep_;

        //- Current number of steps written
        label writeStepI_;

        //- Action to take when conditions are satisfied
        satisfiedAction satisfiedAction_;

        //- Trigger index if satisfiedAction is setTrigger
        label triggerIndex_;

        //- Active flag
        //  Used in the trigger case to bypass any evaluations after the
        //  trigger has been set
        bool active_;

        // Whethere to trigger function objects on stopping
        Switch functionsOnStopping_;

        //- Solution interval - only solve every n outer iterations
        //  Used in the modifiedSolver case
        autoPtr<Function1<scalar>> solveInterval_;

        const word logicalOperator_;

        bool postProcessMode_;

        static const NamedEnum<logicalOps, 2> conditionLogicalOps_;

    // Private Member Functions

        //- No copy construct
        runTimeControl(const runTimeControl&) = delete;

        //- No copy assignment
        void operator=(const runTimeControl&) = delete;


public:

    //- Runtime type information
    TypeName("runTimeControl");


    // Constructors

        //- Construct for given objectRegistry and dictionary
        runTimeControl
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~runTimeControl() = default;


    // Member Functions

        //- Read the runTimeControl data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Calculate the runTimeControl and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimeControls
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
