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
    (c) ICE Stroemungsfoschungs GmbH

Class
    timeDependentFvSolutionFvSchemesFunctionObject

Description
    Set the fvSolution and fvSchemes depending on the simulation time

    The entries to be set are in sub-dictionaries of the original files

    The special name 'reset' resets to the original state
SourceFiles
    timeDependentFvSolutionFvSchemesFunctionObject.C

Contributors/Copyright:
    2015-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef timeDependentFvSolutionFvSchemesFunctionObject_H
#define timeDependentFvSolutionFvSchemesFunctionObject_H

#include "manipulateFvSolutionFvSchemes/manipulateFvSolutionFvSchemesFunctionObject/manipulateFvSolutionFvSchemesFunctionObject.H"
#include "primitives/Tuple2/Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class timeDependentFvSolutionFvSchemesFunctionObject Declaration
\*---------------------------------------------------------------------------*/


class timeDependentFvSolutionFvSchemesFunctionObject
:
    public manipulateFvSolutionFvSchemesFunctionObject
{
    // Private Data

        typedef Tuple2<scalar,word> TriggerTimeName;

        List<TriggerTimeName> solutionTriggers_;
        label currentSolutionTrigger_;

        List<TriggerTimeName> schemesTriggers_;
        label currentSchemesTrigger_;

        bool resetBeforeTrigger_;


    // Private Member Functions

        void checkTriggerList
        (
            const List<TriggerTimeName>& triggers,
            const dictionary& dict
        );

        bool triggerTrigger
        (
            const Time& t,
            const List<TriggerTimeName>& triggers,
            dictionary& dict,
            label& current
        );


public:

    //- Runtime type information
    TypeName("timeDependentFvSolutionFvSchemes");


    // Constructors

        //- Construct from components
        timeDependentFvSolutionFvSchemesFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        timeDependentFvSolutionFvSchemesFunctionObject
        (
            const timeDependentFvSolutionFvSchemesFunctionObject&
        ) = delete;


    // Member Functions

        bool manipulateFvSolution(const Time& t);

        bool manipulateFvSchemes(const Time& t);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const timeDependentFvSolutionFvSchemesFunctionObject&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
