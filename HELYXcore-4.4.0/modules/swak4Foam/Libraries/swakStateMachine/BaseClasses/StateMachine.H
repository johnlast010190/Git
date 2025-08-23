/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : dev
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
    Foam::StateMachine

Description
    Basic statemachine. Stores the state and is in charge of calculating
    progress

    Also holds the ValueExpressionDriver that does the state calculations

SourceFiles
    StateMachine.C

Contributors/Copyright:
    2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef StateMachine_H
#define StateMachine_H

#include "CommonValueExpressionDriver.H"
#include "speciesTable/speciesTable.H"
#include "containers/Lists/PtrList/PtrList.H"
#include "namedEnums/LogicalAccumulationNamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    class StateMachineRepository;

/*---------------------------------------------------------------------------*\
                        Class StateMachine Declaration
\*---------------------------------------------------------------------------*/

class StateMachine
{
    // Private Member Functions

        //- Construct null
        StateMachine(const IOobject&);

        //- Constructor from components
        StateMachine(const dictionary& dict, const fvMesh& mesh);


    // Private data

        //- The driver
        autoPtr<CommonValueExpressionDriver> driver_;

        //- The mesh
        const fvMesh& mesh_;

        //- The state names
        const speciesTable names_;

        //- The name of the machine
        word machineName_;

        //- Initial state
        label initialState_;

        //- Current state
        label state_;

        //- Last time the state changed
        scalar lastStateChange_;

        //- Steps since last change
        label stepsSinceChange_;

        //- How often changed the machine to this state
        labelList changedTo_;

    class StateTransition
    {

        //- State this goes from
        label from_;

        //- State this goes to
        label to_;

        //- The "parent"
        StateMachine &machine_;

        //- Description of the transition
        string description_;

        //- The condition expression
        const exprString condition_;

        //- How should the condition be aqccumulated
        LogicalAccumulationNamedEnum::value logicalAccumulation_;


    public:

        StateTransition
        (
            StateMachine& machine,
            const dictionary& data
        );

        inline label from()
        {
            return from_;
        }

        inline label to()
        {
            return to_;
        }

        inline const string& description()
        {
            return description_;
        }

        bool operator()();
    };

    PtrList<StateTransition> transitions_;


protected:

    // Protected Member Functions

        void resetState
        (
            const word& state,
            scalar timeOfChange,
            label stepsSinceChange,
            labelList changedTo
        );

        inline scalar lastStateChange() const
        {
            return lastStateChange_;
        }


public:

    //- Runtime type information
    TypeName("StateMachine");


    //- Destructor
    virtual ~StateMachine() {}

    //- Construct a new machine and register it
    static StateMachine& NewMachine
    (
        const dictionary& dict,
        const fvMesh& mesh
    );

    static StateMachine& machine
    (
        const word& name
    );

    //- Make sure that the state is consistently written
    static void ensureWrite();

    inline CommonValueExpressionDriver& driver()
    {
        return driver_();
    }

    inline label numberOfStates() const
    {
        return names_.size();
    }

    inline const word& name()
    {
        return machineName_;
    }

    inline label stateCode(word name) const
    {
        return names_[name];
    }

    inline bool hasState(word name) const
    {
        return names_.found(name);
    }

    inline const word &stateName(label stateCode) const
    {
        return names_[stateCode];
    }

    inline label currentState() const
    {
        return state_;
    }

    inline label stepsSinceChange() const
    {
        return stepsSinceChange_;
    }

    inline label changedTo(label state) const
    {
        return changedTo_[state];
    }

    //- Time since last state change
    scalar timeSinceChange() const;

    //- Return value is a string with the information about what happened
    std::string step();

    //- Return value is a string with the information about what happened
    std::string force(label state);

    friend class StateMachineRepository;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const StateMachine&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
