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
    (c) 2024 Engys Ltd.

Class
    setDeltaTBySwakExpressionFunctionObject

Description
    Set the time-step according to a swakExpression

SourceFiles
    setDeltaTBySwakExpressionFunctionObject.C

Contributors/Copyright:
    2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef setDeltaTBySwakExpressionFunctionObject_H
#define setDeltaTBySwakExpressionFunctionObject_H

#include "timeManipulation/timeManipulationFunctionObject/timeManipulationFunctionObject.H"

#include "swakDataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
#ifdef FOAM_DATAENTRY_IS_NOW_FUNCTION1
    using namespace Function1Types;
#endif

/*---------------------------------------------------------------------------*\
                           Class setDeltaTBySwakExpressionFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class setDeltaTBySwakExpressionFunctionObject
:
    public timeManipulationFunctionObject
{
    // Private Data

        swakDataEntry<scalar> deltaTExpression_;


public:

    //- Runtime type information
    TypeName("setDeltaTBySwakExpression");


    // Constructors

        //- Construct from components
        setDeltaTBySwakExpressionFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        setDeltaTBySwakExpressionFunctionObject
        (
            const setDeltaTBySwakExpressionFunctionObject&
        ) = delete;


    // Member Functions

        virtual scalar deltaT();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const setDeltaTBySwakExpressionFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
