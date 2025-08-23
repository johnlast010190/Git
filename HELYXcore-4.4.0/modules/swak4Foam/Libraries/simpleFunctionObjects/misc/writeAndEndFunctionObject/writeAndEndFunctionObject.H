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
    writeAndEndFunctionObject

Description
    General function object for ending a run according to a condition

SourceFiles
    writeAndEndFunctionObject.C

Contributors/Copyright:
    2012-2013, 2015-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef writeAndEndFunctionObject_H
#define writeAndEndFunctionObject_H

#include "simpleFunctionObject/simpleFunctionObject.H"
#include "DataStructures/TimeClone/TimeCloneList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class writeAndEndFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class writeAndEndFunctionObject
:
    public simpleFunctionObject
{
    // Private Data

        //- Flag if the run is already stopped
        bool isStopped_;

        autoPtr<TimeCloneList> lastTimes_;

        bool storeAndWritePreviousState_;


protected:

    // Protected Member Functions

        //- Do the actual checking
        virtual bool endRunNow() = 0;

        bool isStopped() const
        {
            return isStopped_;
        }


public:

    //- Runtime type information
    TypeName("writeAndEnd");


    // Constructors

        //- Construct from components
        writeAndEndFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        writeAndEndFunctionObject(const writeAndEndFunctionObject&) = delete;


    // Member Functions

        bool start();

        void writeSimple();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const writeAndEndFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
