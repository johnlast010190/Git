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
    Foam::swakCodedFunctionObject

Description
    functionObject using dynamic code compilation. Reads and writes swak
    global variables

SourceFiles
    swakCodedFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef swakCodedFunctionObject_H
#define swakCodedFunctionObject_H

#include "codedFunctionObject/codedFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                    Class swakCodedFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class swakCodedFunctionObject
:
    public codedFunctionObject
{
    // Private Data

        //- Variables in those namespaces should be inserted into the
        //  Coded namespace before running the snipplet
        wordList swakToCodedNamespaces_;

        //- Namespace to write values to
        //  This global swak-namespace will get all values (if the user request such)
        word codedToSwakNamespace_;

        //- Coded variables that swak will receive
        //  The values of these variables will get written to the global swak namespace
        //  only variables that can be either interpreted as a scalar or a vector are
        //  supported
        wordList codedToSwakVariables_;

        //- Should the code have Debug-Information
        bool verboseCode_;


protected:

    //- Insert swak-specific code
    void injectSwakCode(const word &key,const string &pre,const string &post);


public:

    //- Runtime type information
    TypeName("swakCoded");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        swakCodedFunctionObject
        (
            const word& name,
            const Time& time,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        swakCodedFunctionObject(const swakCodedFunctionObject&) = delete;


    //- Destructor
    virtual ~swakCodedFunctionObject();


    // Member Functions

        //- Read and set the function object if its data have changed
        virtual bool read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const swakCodedFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
