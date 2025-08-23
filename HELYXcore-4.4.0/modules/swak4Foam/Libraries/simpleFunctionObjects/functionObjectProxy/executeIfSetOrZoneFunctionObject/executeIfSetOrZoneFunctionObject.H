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
    executeIfSetOrZoneFunctionObject

Description
    This proxy-functionObject only executes its childrens if the sets and zones
    listed in a dictionary exist

SourceFiles
    executeIfSetOrZoneFunctionObject.C

Contributors/Copyright:
    2011, 2013-2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef executeIfSetOrZoneFunctionObject_H
#define executeIfSetOrZoneFunctionObject_H

#include "functionObjectProxy/conditionalFunctionObjectListProxy/conditionalFunctionObjectListProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
            Class executeIfSetOrZoneFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class executeIfSetOrZoneFunctionObject
:
    public conditionalFunctionObjectListProxy
{
    // Private Data

        //- The dictionary with the sets and zones
        dictionary setsAndZones_;

        const polyMesh& mesh_;

        bool loadAndCacheMissingSets_;


    // Private Member Functions

        //- Reread the regular expression
        void readSetsAndZones(const dictionary& dict);

        //- Check if set is there and load if missing
        template<class T>
        bool hasSet(const word& name);


protected:

    // Protected Member Functions

        //- The actual conition that determines
        //  whether the FO should be executed
        virtual bool condition();


public:

    //- Runtime type information
    TypeName("executeIfSetOrZone");


    // Constructors

        //- Construct from components
        executeIfSetOrZoneFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        executeIfSetOrZoneFunctionObject
        (
            const executeIfSetOrZoneFunctionObject&
        ) = delete;


    // Member Functions

        bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const executeIfSetOrZoneFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
