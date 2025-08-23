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
    stateMachineStateFunctionObject

Description
    Writes the state of a state machine into a file

SourceFiles
    stateMachineStateFunctionObject.C

Contributors/Copyright:
    2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#ifndef stateMachineStateFunctionObject_H
#define stateMachineStateFunctionObject_H

#include "timelineFunctionObject/timelineFunctionObject.H"

#include "BaseClasses/StateMachine.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class stateMachineStateFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class stateMachineStateFunctionObject
:
    public timelineFunctionObject
{
    //- Evaluate the expression
    template <class T>
    Field<T> evaluate(const word& fieldName,T preset) const;

    //- Name of the state machinr
    word machineName_;


protected:

    word dirName();

    stringList columnNames();

    wordList fileNames();


public:

    //- Runtime type information
    TypeName("stateMachineState");


    // Constructors

        //- Construct from components
        stateMachineStateFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        stateMachineStateFunctionObject
        (
            const stateMachineStateFunctionObject&
        ) = delete;


    //- Do the writing
    virtual void writeSimple();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const stateMachineStateFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
