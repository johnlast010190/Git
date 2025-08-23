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
    stateMachineCreateAndUpdateFunctionObject

Description
    Creates a state-machine and updates it at ever time-step

SourceFiles
    stateMachineCreateAndUpdateFunctionObject.C

Contributors/Copyright:
    2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef stateMachineCreateAndUpdateFunctionObject_H
#define stateMachineCreateAndUpdateFunctionObject_H

#include "simpleFunctionObject/simpleFunctionObject.H"

#include "BaseClasses/StateMachine.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class stateMachineCreateAndUpdateFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class stateMachineCreateAndUpdateFunctionObject
:
    public simpleFunctionObject
{
    // Private Data

        StateMachine& machine_;


public:

    //- Runtime type information
    TypeName("stateMachineCreateAndUpdate");


    // Constructors

        //- Construct from components
        stateMachineCreateAndUpdateFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        stateMachineCreateAndUpdateFunctionObject
        (
            const stateMachineCreateAndUpdateFunctionObject&
        ) = delete;


    // Member Functions

        bool start();

        void writeSimple();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const stateMachineCreateAndUpdateFunctionObject&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
