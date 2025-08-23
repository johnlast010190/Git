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

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "ConcretePluginFunction.H"

namespace Foam {


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class DriverType>
ConcretePluginFunction<DriverType>::ConcretePluginFunction(
    const DriverType &parentDriver,
    const word &name,
    const word &returnType,
    const string &argumentSpecification
):
    CommonPluginFunction(
        parentDriver,
        name,
        returnType,
        argumentSpecification
    )
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class DriverType>
autoPtr<ConcretePluginFunction<DriverType>> ConcretePluginFunction<DriverType>::New (
    const DriverType& driver,
    const word &name
)
{
    if (debug) {
        Info<< "ConcretePluginFunction::New looking for "
            << name
            << endl;
    }
    typename nameConstructorTable::iterator cstrIter =
        nameConstructorTable_().find(name);
    if (cstrIter==nameConstructorTable_().end()) {
        FatalErrorIn(
            "autoPtr<ConcretePluginFunction> ConcretePluginFunction<"+
            DriverType::typeName+">::New"
        ) << "Unknown plugin function " << name << endl
                << exit(FatalError);
    }
    return autoPtr<ConcretePluginFunction>
        (
            cstrIter->second(driver,name)
        );
}

template<class DriverType>
bool ConcretePluginFunction<DriverType>::exists (
    const DriverType& driver,
    const word &name
)
{
    static bool firstCall=true;
    if (firstCall) {
        firstCall=false;

        if (nameConstructorTable_().size()>0) {
            Info<< endl << "Loaded plugin functions for '"+
                DriverType::typeName+"':" << endl;

            for (const auto& iter : nameConstructorTable_()) {
                const word &theName = iter.first;
                Info<< "  " << theName << ":" << endl
                    << "    " << iter.second(driver,theName)->helpText() << endl;
            }

            Info<< endl;
        }
    }

    typename nameConstructorTable::iterator cstrIter =
        nameConstructorTable_().find(name);

    return cstrIter!=nameConstructorTable_().end();
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace

// ************************************************************************* //
