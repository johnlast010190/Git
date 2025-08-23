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
    simpleFunctionObject

Description
    Basis for the other function objects here. Creates the directory and the files

SourceFiles
    simpleFunctionObject.C

Contributors/Copyright:
    2008-2011, 2013, 2015-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef simpleDataFunctionObject_H
#define simpleDataFunctionObject_H

#include "simpleFunctionObject/simpleFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class simpleDataFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class simpleDataFunctionObject
:
    public simpleFunctionObject
{
    // Private Data

        static fileName defaultPostProcDir_;

        fileName postProcDir_;


protected:

    // Protected Member Functions

        //- Get the path to the data directory
        virtual fileName dataDir();

        //- Get the path to the base directory
        fileName baseDir();

        //- Name of the sub-dictionary of the case
        virtual word dirName()=0;


public:

    //- Runtime type information
    TypeName("simpleDataFunctionObject");


    // Constructors

        //- Construct from components
        simpleDataFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        simpleDataFunctionObject(const simpleDataFunctionObject&) = delete;


    // Member Functions

        bool start();

        static void setPostProcDir(const fileName &);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const simpleDataFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
