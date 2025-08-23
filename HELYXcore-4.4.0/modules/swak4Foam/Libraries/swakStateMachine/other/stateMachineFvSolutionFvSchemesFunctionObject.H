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
    stateMachineFvSolutionFvSchemesFunctionObject

Description
    Set the fvSolution and fvSchemes depending on the state of a state machine

    The entries to be set are in sub-dictionaries of the original files

    States that are not in the mapping mean that the default values are used

SourceFiles
    stateMachineFvSolutionFvSchemesFunctionObject.C

Contributors/Copyright:
    2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef stateMachineFvSolutionFvSchemesFunctionObject_H
#define stateMachineFvSolutionFvSchemesFunctionObject_H

#include "manipulateFvSolutionFvSchemes/manipulateFvSolutionFvSchemesFunctionObject/manipulateFvSolutionFvSchemesFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class stateMachineFvSolutionFvSchemesFunctionObject Declaration
\*---------------------------------------------------------------------------*/


class stateMachineFvSolutionFvSchemesFunctionObject
:
    public manipulateFvSolutionFvSchemesFunctionObject
{
    // Private Data

        typedef HashTable<word,word> StateToMap;

        const word solutionStateMachineName_;
        const StateToMap stateToSolutionMapping_;
        word lastSolutionState_;

        const word schemesStateMachineName_;
        const StateToMap stateToSchemesMapping_;
        word lastSchemesState_;

        bool resetBeforeTrigger_;


    // Private Member Functions

        void checkTriggerMapping
        (
            const word& machine,
            const StateToMap& mapping,
            const dictionary& dict
        );

        bool triggerMachine
        (
            const Time& t,
            const word& machine,
            const StateToMap& mapping,
            dictionary& dict,
            word& current
        );


public:

    //- Runtime type information
    TypeName("stateMachineFvSolutionFvSchemes");


    // Constructors

        //- Construct from components
        stateMachineFvSolutionFvSchemesFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        stateMachineFvSolutionFvSchemesFunctionObject
        (
            const stateMachineFvSolutionFvSchemesFunctionObject&
        ) = delete;


    // Member Functions

        bool manipulateFvSolution(const Time& t);

        bool manipulateFvSchemes(const Time& t);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const stateMachineFvSolutionFvSchemesFunctionObject&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
