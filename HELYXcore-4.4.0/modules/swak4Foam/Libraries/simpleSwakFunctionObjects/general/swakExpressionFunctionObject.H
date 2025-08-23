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
    swakExpressionFunctionObject

Description
    Can evaluate any swak-expression that can be produced by runTime-selection

SourceFiles
    swakExpressionFunctionObject.C

Contributors/Copyright:
    2010-2011, 2013-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#ifndef swakExpressionFunctionObject_H
#define swakExpressionFunctionObject_H

#include "timelineFunctionObject/timelineFunctionObject.H"

#include "CommonValueExpressionDriver.H"

#include "namedEnums/NumericAccumulationNamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class swakExpressionFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class swakExpressionFunctionObject
:
    public timelineFunctionObject
{
    // Private Data

        //- The expression that is to be evaluated
        exprString expression_;

        //- The mask expression
        exprString maskExpression_;

        //- List with the accumulations of evaluated expression
        List<NumericAccumulationNamedEnum::accuSpecification> accumulations_;

        //- The expression driver
        autoPtr<CommonValueExpressionDriver> driver_;


    // Private Member Functions

        //- Evaluate the expression
        template <class T>
        Field<T> evaluate(const word& fieldName, T preset) const;


protected:

    // Protected Member Functions

        word dirName();

        stringList columnNames();

        wordList fileNames();

        template<class T>
        void writeTheData(CommonValueExpressionDriver &driver);


public:

    //- Runtime type information
    TypeName("swakExpression");


    // Constructors

        //- Construct from components
        swakExpressionFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        swakExpressionFunctionObject
        (
            const swakExpressionFunctionObject&
        ) = delete;


    // Member Functions

        //- Do the writing
        virtual void writeSimple();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const swakExpressionFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "swakExpressionFunctionObjectTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
