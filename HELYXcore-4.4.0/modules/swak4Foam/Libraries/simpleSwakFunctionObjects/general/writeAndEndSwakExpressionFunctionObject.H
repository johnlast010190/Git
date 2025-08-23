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
    writeAndEndSwakExpressionFunctionObject

Description
    Check whether a swak-condition is met. If yes everything is written to disk
    and the run ends gracefully

SourceFiles
    writeAndEndSwakExpressionFunctionObject.C

Contributors/Copyright:
    2012-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#ifndef writeAndEndSwakExpressionFunctionObject_H
#define writeAndEndSwakExpressionFunctionObject_H

#include "misc/writeAndEndFunctionObject/writeAndEndFunctionObject.H"

#include "CommonValueExpressionDriver.H"

#include "namedEnums/LogicalAccumulationNamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class writeAndEndSwakExpressionFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class writeAndEndSwakExpressionFunctionObject
:
    public writeAndEndFunctionObject
{
    // Private Data

        //- The expression driver
        autoPtr<CommonValueExpressionDriver> driver_;

        //- The logical expression
        exprString logicalExpression_;

        //- Possible logical accumulations
        LogicalAccumulationNamedEnum::value logicalAccumulation_;


    // Private Member Functions

        void readParameters(const dictionary &);


protected:

    // Protected Member Functions

        //- Do the actual checking
        virtual bool endRunNow();


public:

    //- Runtime type information
    TypeName("writeAndEndSwakExpression");


    // Constructors

        //- Construct from components
        writeAndEndSwakExpressionFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        writeAndEndSwakExpressionFunctionObject
        (
            const writeAndEndSwakExpressionFunctionObject&
        ) = delete;


    // Member Functions

        bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const writeAndEndSwakExpressionFunctionObject&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
