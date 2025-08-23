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
    executeIfEnvironmentVariableFunctionObject

Description
    This proxy-functionObject only executs its childrens if an environment
    variable fits a certain condition:
    - exists
    - does not exist
    - fits a regular expression

SourceFiles
    executeIfEnvironmentVariableFunctionObject.C

Contributors/Copyright:
    2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef executeIfEnvironmentVariableFunctionObject_H
#define executeIfEnvironmentVariableFunctionObject_H

#include "functionObjectProxy/conditionalFunctionObjectListProxy/conditionalFunctionObjectListProxy.H"
#include "regExp.H"
#include "primitives/enums/NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class executeIfEnvironmentVariableFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class executeIfEnvironmentVariableFunctionObject
:
    public conditionalFunctionObjectListProxy
{
    // Private Data

        //- Name of the variable
        word variableName_;

        //- The regular expression against which the content
        //  of the environment variable should be matched
        regExp contentRegexp_;


    // Private Member Functions

        //- Reread the data
        void readData(const dictionary& dict);


public:

    // some compilers need that
    enum fitVariableMode
    {
        fvmExists,
        fvmDoesNotExist,
        fvmFitsRegexp
    };


private:

    // Private Data

        static const NamedEnum<fitVariableMode, 3> fitVariableModeNames_;

        //- The mode with which the variable should be fitted
        fitVariableMode fitVariableMode_;


protected:

    // Protected Member Functions

        //- The actual conition that determines
        //  whether the FO should be executed
        virtual bool condition();


public:

    //- Runtime type information
    TypeName("executeIfEnvironmentVariable");


    // Constructors

        //- Construct from components
        executeIfEnvironmentVariableFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        executeIfEnvironmentVariableFunctionObject
        (
            const executeIfEnvironmentVariableFunctionObject&
        ) = delete;


    // Member Functions

        bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const executeIfEnvironmentVariableFunctionObject&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
