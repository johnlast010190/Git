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
    provokeSignalFunctionObject

Description
    Raises a specified signal at a specified time. For parallel runs also
    the processors on which to raise the signal can be specified.

    The only purpose of this function object is the testing of the error
    handling

SourceFiles
    provokeSignalFunctionObject.C

Contributors/Copyright:
    2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef provokeSignalFunctionObject_H
#define provokeSignalFunctionObject_H

#include "simpleFunctionObject/simpleFunctionObject.H"
#include "primitives/enums/NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class provokeSignalFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class provokeSignalFunctionObject
:
    public simpleFunctionObject
{
public:

    enum possibleSignals
    {
        sigFPE,
        sigSEGV,
        sigINT,
        sigTERM,
        sigQUIT,
        sigUSR1,
        sigUSR2
    };


private:

    // Private Data

        possibleSignals signalToRaise_;

        scalar timeToRaise_;

        bool raiseOnThisProc_;

        static const NamedEnum<possibleSignals, 7> possibleSignalsNames_;


public:

    //- Runtime type information
    TypeName("provokeSignal");


    // Constructors

        //- Construct from components
        provokeSignalFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        provokeSignalFunctionObject
        (
            const provokeSignalFunctionObject&
        ) = delete;


    // Member Functions

        void writeSimple();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const provokeSignalFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
