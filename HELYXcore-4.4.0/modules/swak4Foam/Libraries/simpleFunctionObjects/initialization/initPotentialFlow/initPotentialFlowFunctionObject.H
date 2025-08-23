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
    initPotentialFlowFunctionObject

Description
    Calculate the potential flow field solution for a field. Either just at
    the start or at every time step

SourceFiles
    initPotentialFlowFunctionObject.C

Contributors/Copyright:
    2008-2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef initPotentialFlowFunctionObject_H
#define initPotentialFlowFunctionObject_H

#include "update/updateSimple/updateSimpleFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class initPotentialFlowFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class initPotentialFlowFunctionObject
:
    public updateSimpleFunctionObject
{
    // Private Data

        //- Name of the velocity field
        word UName_;

        //- Name of the pressure
        word pName_;

        //- Write the fields?
        bool writeFields_;

        //- Write the old fields?
        bool writeOldFields_;

        //- Override pressure?
        bool overrideP_;

        //- Initialize the U BCs?
        bool initialiseUBCs_;


protected:

    //- do the actual potential flow calculation
    void recalc();


public:

    //- Runtime type information
    TypeName("initPotentialFlow");

    // Constructors

        //- Construct from components
        initPotentialFlowFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        initPotentialFlowFunctionObject
        (
            const initPotentialFlowFunctionObject&
        ) = delete;


    // Member Functions

        bool start();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const initPotentialFlowFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
