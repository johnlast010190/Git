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
    (c) 1991-2008 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::CommonPluginFunction

Description

SourceFiles
    CommonPluginFunctionI.H
    CommonPluginFunction.C
    CommonPluginFunctionIO.C

Contributors/Copyright:
    2012-2013, 2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef CommonPluginFunction_H
#define CommonPluginFunction_H

#include "primitives/strings/lists/wordList.H"
#include "primitives/strings/lists/stringList.H"

#include "db/typeInfo/typeInfo.H"

#include "ExpressionResult/ExpressionResult.H"

#include "primitives/Vector/vector/vector.H"
#include "primitives/Tensor/tensor/tensor.H"
#include "primitives/SymmTensor/symmTensor/symmTensor.H"
#include "primitives/SphericalTensor/sphericalTensor/sphericalTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    class CommonValueExpressionDriver;
    class fvMesh;

/*---------------------------------------------------------------------------*\
                    Class CommonPluginFunction Declaration
\*---------------------------------------------------------------------------*/
class CommonPluginFunction
{
    // Private Data

        //- The driver that calls this plugin
        const CommonValueExpressionDriver &parentDriver_;

        //- Name of the function
        word name_;

        //- Type of data it returns (textual representation)
        word returnType_;

        //- Names of the arguments
        wordList argumentNames_;

        //- Which parser can read the arguments
        wordList argumentParsers_;

        //- Expected types of the arguments
        wordList argumentTypes_;

        //- The actual expression result. To be set by doEvaluation
        ExpressionResult result_;


    // Private Member Functions

        void evaluateInternal(const string& content, int& consumed);


protected:

    // Protected Member Functions

        //- Wrapper for argument reading
        label readArgument
        (
            label index,
            const string& content,
            const word& type
        );

        //- Set a primitive argument
        virtual void setArgument(label index, const string& value)
        {
            // should be overwritten by the function itself
            notImplemented("setArgument(label index,const &string)");
        }

        //- Set a primitive argument
        virtual void setArgument(label index, const word& value)
        {
            // should be overwritten by the function itself
            notImplemented("setArgument(label index,const &word)");
        }

        //- Set a primitive argument
        virtual void setArgument(label index, const scalar& value)
        {
            // should be overwritten by the function itself
            notImplemented("setArgument(label index,const &scalar)");
        }

        //- Set a primitive argument
        virtual void setArgument(label index, const label& value)
        {
            // should be overwritten by the function itself
            notImplemented("setArgument(label index,const &label)");
        }

        //- Set a primitive argument
        virtual void setArgument(label index, const vector& value)
        {
            // should be overwritten by the function itself
            notImplemented("setArgument(label index,const &vector)");
        }

        //- set a primitive argument
        virtual void setArgument(label index, const tensor& value)
        {
            // should be overwritten by the function itself
            notImplemented("setArgument(label index,const &tensor)");
        }

        //- Set a primitive argument
        virtual void setArgument(label index, const symmTensor& value)
        {
            // should be overwritten by the function itself
            notImplemented("setArgument(label index,const &symmTensor)");
        }

        //- Set a primitive argument
        virtual void setArgument(label index, const sphericalTensor& value)
        {
            // should be overwritten by the function itself
            notImplemented("setArgument(label index,const &sphericalTensor)");
        }

        //- Set a primitive argument
        virtual void setArgument(label index, const bool& value)
        {
            // should be overwritten by the function itself
            notImplemented("setArgument(label index,const &bool)");
        }

        //- Wrapper for argument reading
        label readArgument
        (
            label index,
            const string& content,
            CommonValueExpressionDriver& driver
        );

        //- Read an argument that is to be parsed
        virtual void setArgument
        (
            label index,
            const string& content,
            const CommonValueExpressionDriver& driver
        )
        {
            // should be overwritten by the function itself
            notImplemented("setArgument(label index,const string &content,const CommonValueExpressionDriver&)");
        }

        //- Get an addtional ID for the parser of argument i
        virtual word getID(label index)
        {
            return "method_getID_notImplemented";
        }

        //- Do the actual evaluation
        virtual void doEvaluation() = 0;

        //- Scan if no values are bein read
        label scanEmpty(const string& content, word sym);

        inline CommonValueExpressionDriver& parentDriver();

        inline const fvMesh& mesh() const;

        inline ExpressionResult& result();


public:

    //- Runtime type information
    TypeName("CommonPluginFunction");


    // Constructors

        //- Construct from components
        CommonPluginFunction
        (
            const CommonValueExpressionDriver& parentDriver,
            const word& name,
            const word& returnType,
            const string& argumentSpecificationString
        );

        //- Disallow default bitwise copy construction
        CommonPluginFunction(const CommonPluginFunction&) = delete;


    //- Destructor
    virtual ~CommonPluginFunction()
    {}



    // Public Member Functions

        //- Return the name
        const word& name() const
        {
            return name_;
        }

        const word& returnType() const
        {
            return returnType_;
        }

        word& returnType()
        {
            return returnType_;
        }

        //- Return a help text constructed from parameters
        string helpText() const;

        //- Evaluate the string and report how many characters were consumed in total
        template <typename T>
        tmp<T> evaluate(const string& content, int& consumed);

        template <typename T>
        tmp<Field<T>> evaluate
        (
            const string& content,
            int& consumed,
            bool isPoint
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const CommonPluginFunction&) = delete;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CommonPluginFunctionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
