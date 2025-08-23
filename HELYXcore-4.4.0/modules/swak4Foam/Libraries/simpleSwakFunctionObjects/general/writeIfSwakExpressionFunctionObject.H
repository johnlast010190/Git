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
    writeIfSwakExpressionFunctionObject

Description

    Check whether a swak-condition is met. If yes start writing
    again. Also offers the possibility to stop writing via an
    expression

SourceFiles
    writeIfSwakExpressionFunctionObject.C

Contributors/Copyright:
    2012-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#ifndef writeIfSwakExpressionFunctionObject_H
#define writeIfSwakExpressionFunctionObject_H

#include "writeFields/conditionDrivenWriting/conditionDrivenWritingFunctionObject.H"

#include "CommonValueExpressionDriver.H"

#include "namedEnums/LogicalAccumulationNamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class writeIfSwakExpressionFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class writeIfSwakExpressionFunctionObject
:
    public conditionDrivenWritingFunctionObject
{
    // Private Data

        //- The expression driver
        autoPtr<CommonValueExpressionDriver> driver_;

        //- The logical expression
        exprString writeExpression_;

        //- Possible logical accumulations
        LogicalAccumulationNamedEnum::value writeAccumulation_;

        //- The logical expression
        exprString stopWriteExpression_;

        //- Possible logical accumulations
        LogicalAccumulationNamedEnum::value stopWriteAccumulation_;

        //- The logical expression
        exprString stopCooldownExpression_;

        //- Possible logical accumulations
        LogicalAccumulationNamedEnum::value stopCooldownAccumulation_;


    // Private Member Functions

        void readParameters(const dictionary&);

        bool evaluateCondition
        (
            exprString expression,
            LogicalAccumulationNamedEnum::value accumulation
        );


protected:

    // Protected Member Functions

        //- Should writing start now?
        virtual bool checkStartWriting();

        //- Switch off writing?
        virtual bool checkStopWriting();

        //- Switch off cooldown?
        virtual bool checkStopCooldown();


public:

    //- Runtime type information
    TypeName("writeIfSwakExpression");


    // Constructors

        //- Construct from components
        writeIfSwakExpressionFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        writeIfSwakExpressionFunctionObject
        (
            const writeIfSwakExpressionFunctionObject&
        ) = delete;


    // Member Functions

        bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const writeIfSwakExpressionFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
