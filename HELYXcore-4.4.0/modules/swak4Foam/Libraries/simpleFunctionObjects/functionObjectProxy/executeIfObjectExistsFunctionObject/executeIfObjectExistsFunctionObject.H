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
    executeIfObjectExistsFunctionObject

Description
    This proxy-functionObject only executs its childrens if a certain
    object exists (or does not exist) in the registry

SourceFiles
    executeIfObjectExistsFunctionObject.C

Contributors/Copyright:
    2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef executeIfObjectExistsFunctionObject_H
#define executeIfObjectExistsFunctionObject_H

#include "functionObjectProxy/conditionalFunctionObjectListProxy/conditionalFunctionObjectListProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class executeIfObjectExistsFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class executeIfObjectExistsFunctionObject
:
    public conditionalFunctionObjectListProxy
{
    // Private Data

        //- Name of the object to look for
        word objectName_;

        //- Should the type be checked or ist the name sufficient
        bool checkType_;

        //- Type of the object
        word objectType_;

        //- Should the object exist or not
        bool objectShouldExist_;


    // Private Member Functions

        void readParameters(const dictionary &);


protected:

    // Protected Member Functions

        //- The actual conition that determines
        //  whether the FO should be executed
        virtual bool condition();


public:

    //- Runtime type information
    TypeName("executeIfObjectExists");


    // Constructors

        //- Construct from components
        executeIfObjectExistsFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        executeIfObjectExistsFunctionObject
        (
            const executeIfObjectExistsFunctionObject&
        ) = delete;


    // Member Functions

        bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const executeIfObjectExistsFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
