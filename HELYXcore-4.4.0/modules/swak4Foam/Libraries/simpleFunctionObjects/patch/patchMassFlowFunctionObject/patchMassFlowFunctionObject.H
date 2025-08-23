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
    patchMassFlowFunctionObject

Description
    Calculates the mass flow on the patches

SourceFiles
    patchMassFlowFunctionObject.C

Contributors/Copyright:
    2008-2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef patchMassFlowFunctionObject_H
#define patchMassFlowFunctionObject_H

#include "patch/patchFunctionObject/patchFunctionObject.H"
#include "SolverInfo/SolverInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class patchMassFlowFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class patchMassFlowFunctionObject
:
    public patchFunctionObject
{
    // Private Data

        SolverInfo solver_;


protected:

    // Protected Member Functions

        //- Names of the files
        virtual wordList fileNames();

        word dirName();


public:

    //- Runtime type information
    TypeName("patchMassFlow");


    // Constructors

        //- Construct from components
        patchMassFlowFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        patchMassFlowFunctionObject
        (
            const patchMassFlowFunctionObject&
        ) = delete;


    // Member Functions

        //- Sample and write
        void writeSimple();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const patchMassFlowFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
