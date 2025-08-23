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
    (c) held by original author
    (c) 2024 Engys Ltd.

Class
    executeIfSwakExpressionFunctionObject

Description
    This proxy-functionObject only executs its childrens if a
    swak-expression evaluates to true

SourceFiles
    executeIfSwakExpressionFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef executeIfSwakExpressionFunctionObject_H
#define executeIfSwakExpressionFunctionObject_H

#include "functionObjectProxy/conditionalFunctionObjectListProxy/conditionalFunctionObjectListProxy.H"

#include "CommonValueExpressionDriver.H"

#include "namedEnums/LogicalAccumulationNamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class executeIfSwakExpressionFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class executeIfSwakExpressionFunctionObject
:
    public conditionalFunctionObjectListProxy
{
    // Private Data

        //- The expression driver
        autoPtr<CommonValueExpressionDriver> driver_;

        //- The logical expression
        exprString logicalExpression_;

        LogicalAccumulationNamedEnum::value logicalAccumulation_;


    // Private Member Functions

        void readParameters(const dictionary &);


protected:

    // Protected Member Functions

        //- The actual conition that determines
        //  whether the FO should be executed
        virtual bool condition();


public:

    //- Runtime type information
    TypeName("executeIfSwakExpression");


    // Constructors

        //- Construct from components
        executeIfSwakExpressionFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        executeIfSwakExpressionFunctionObject
        (
            const executeIfSwakExpressionFunctionObject&
        ) = delete;


    // Member Functions

        bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const executeIfSwakExpressionFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
