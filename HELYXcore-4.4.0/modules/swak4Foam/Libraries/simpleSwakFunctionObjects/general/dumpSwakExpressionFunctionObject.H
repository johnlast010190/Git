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
    dumpSwakExpressionFunctionObject

Description
    Can evaluate any swak-expression that can be produced by runTime-selection

SourceFiles
    dumpSwakExpressionFunctionObject.C

Contributors/Copyright:
    2011, 2013-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#ifndef dumpSwakExpressionFunctionObject_H
#define dumpSwakExpressionFunctionObject_H

#include "timelineFunctionObject/timelineFunctionObject.H"

#include "CommonValueExpressionDriver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class dumpSwakExpressionFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class dumpSwakExpressionFunctionObject
:
    public timelineFunctionObject
{
    // Private Data

        //- The expression that is to be evaluated
        exprString expression_;

        //- The expression driver
        autoPtr<CommonValueExpressionDriver> driver_;


    // Private Member Functions

        //- Evaluate the expression
        template <class T>
        Field<T> evaluate(const word& fieldName,T preset) const;

        template<class Type>
        static void writeValue(Ostream& o, const Type& val, unsigned int& w);


protected:

    // Protected Member Functions

        word dirName();

        stringList columnNames();

        wordList fileNames();

        template<class T>
        void writeTheData(CommonValueExpressionDriver& driver);


public:

    //- Runtime type information
    TypeName("dumpSwakExpression");


    // Constructors

        //- Construct from components
        dumpSwakExpressionFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        dumpSwakExpressionFunctionObject
        (
            const dumpSwakExpressionFunctionObject&
        ) = delete;


    // Member Functions

        //- Do the writing
        virtual void writeSimple();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const dumpSwakExpressionFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "dumpSwakExpressionFunctionObjectTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
