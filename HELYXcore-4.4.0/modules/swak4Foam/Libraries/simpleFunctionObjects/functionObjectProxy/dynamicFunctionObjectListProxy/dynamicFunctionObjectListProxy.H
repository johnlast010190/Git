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
    dynamicFunctionObjectListProxy

Description
    Dynamically generates a functionObjectList and executes it

SourceFiles
    dynamicFunctionObjectListProxy.C

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef dynamicFunctionObjectListProxy_H
#define dynamicFunctionObjectListProxy_H

#include "functionObjectProxy/functionObjectListProxy/functionObjectListProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class dynamicFunctionObjectListProxy Declaration
\*---------------------------------------------------------------------------*/

class dynamicFunctionObjectListProxy
:
    public functionObjectListProxy
{
public:

    class dynamicDictionaryProvider
    {
    public:

        //- Runtime type information
        TypeName("dynamicDictionaryProvider");

        declareRunTimeSelectionTable
        (
            autoPtr,
            dynamicDictionaryProvider,
            dictionary,
            (
                const dictionary& dict,
                const dynamicFunctionObjectListProxy &owner
            ),
            (dict,owner)
        );


        // Constructors

            //- Construct from components
            dynamicDictionaryProvider
            (
                const dictionary&,
                const dynamicFunctionObjectListProxy &
            )
            {}

            //- Disallow default bitwise copy construction
            dynamicDictionaryProvider
            (
                const dynamicDictionaryProvider&
            ) = delete;

            static autoPtr<dynamicDictionaryProvider> New
            (
                const word& type,
                const dictionary& dict,
                const dynamicFunctionObjectListProxy &owner
            );


        //- Destructor
        virtual ~dynamicDictionaryProvider()
        {}


        // Member Functions

            //- Return the dictionary text
            virtual string getDictionaryText() = 0;


        // Member Operators

            //- Disallow default bitwise assignment
            void operator=(const dynamicDictionaryProvider&) = delete;
    };


protected:

    // Protected Member Functions

        //- initialize the function-objects
        virtual void initFunctions();


    // Protected Data

        //- the generated dictionary
        autoPtr<dictionary> dynamicDict_;

        //- the actual provider of the dictionary text
        autoPtr<dynamicDictionaryProvider> provider_;


public:

    //- Runtime type information
    TypeName("dynamicFunctionObjectListProxy");


    // Constructors

        //- Construct from components
        dynamicFunctionObjectListProxy
        (
            const word&,
            const Time&,
            const dictionary&,
            const char* providerName = ""
        );

        //- Disallow default bitwise copy construction
        dynamicFunctionObjectListProxy
        (
            const dynamicFunctionObjectListProxy&
        ) = delete;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const dynamicFunctionObjectListProxy&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
