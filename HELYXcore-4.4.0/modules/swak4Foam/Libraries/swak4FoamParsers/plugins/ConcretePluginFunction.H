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
    Foam::ConcretePluginFunction

Description

SourceFiles
    ConcretePluginFunction.C

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef ConcretePluginFunction_H
#define ConcretePluginFunction_H

#include "CommonPluginFunction.H"

#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ConcretePluginFunction Declaration
\*---------------------------------------------------------------------------*/
template<class DriverType>
class ConcretePluginFunction
:
    public CommonPluginFunction
{

public:

    //- Runtime type information
    TypeName("ConcretePluginFunction");


    // Constructors

        //- Construct from components
        ConcretePluginFunction
        (
            const DriverType& parentDriver,
            const word& name,
            const word& returnType,
            const string& argumentSpecification
        );

        //- Disallow default bitwise copy construction
        ConcretePluginFunction(const ConcretePluginFunction&) = delete;


    //- Destructor
    virtual ~ConcretePluginFunction() {}


    typedef DriverType driverType;

    declareRunTimeSelectionTable
    (
        autoPtr,
        ConcretePluginFunction,
        name,
        (
            const DriverType& driver,
            const word &name
        ),
        (driver,name)
    );

    // Member Functions

        static bool exists(const DriverType& driver, const word& name);

        static autoPtr<ConcretePluginFunction> New
        (
            const DriverType& driver,
            const word& name
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ConcretePluginFunction&) = delete;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ConcretePluginFunction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
