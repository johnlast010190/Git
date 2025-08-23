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

Contributors/Copyright:
    2011-2013, 2015-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "functionObjectListProxy.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "meshes/polyMesh/polyMesh.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "include/swakTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(functionObjectListProxy, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        functionObjectListProxy,
        dictionary
    );


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

functionObjectListProxy::functionObjectListProxy
(
    const word& name,
    const Time& t,
    const dictionary& dict,
    bool allowReadingDuringConstruction
)
:
    simpleFunctionObject(name, t, dict)
{
    Pbug << "Constructing " << name << endl;

    if
    (
        allowReadingDuringConstruction
     && !dict.found("functions")
    )
    {
        FatalErrorIn("functionObjectListProxy::functionObjectListProxy")
            << "No entry 'functions' in dictionary of " << name << endl
            << exit(FatalError);
    }
    if
    (
        allowReadingDuringConstruction
     && dict.lookup<bool>("readDuringConstruction")
    ) {
        Dbug<< this->name() << " list initialized during construction" << endl;
        read(dict);
    }
    else if
    (
        !allowReadingDuringConstruction
     && dict.lookup<bool>("readDuringConstruction")
    )
    {
        WarningIn("functionObjectListProxy::functionObjectListProxy")
            << "For " << name << " the 'readDuringConstruction'-variable "
            << " is ignored because this concrete class does not allow "
            << " constructing the functionObjects during construction"<< endl;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void functionObjectListProxy::initFunctions()
{
    functions_.set(new functionObjectList(time(), dict_));

    Dbug<< this->name() << " list initialized with "
        << functions_->size() << " FOs" << endl;
}

functionObjectList& functionObjectListProxy::functions()
{
    if (!functions_.valid())
    {
        Dbug<< this->name() << " list initialized on demand" << endl;
        initFunctions();
    }

    return functions_();
}


bool functionObjectListProxy::execute(const bool forceWrite)
{
    Dbug<< this->name() << " functionObjectListProxy::execute()" << endl;

#ifdef FOAM_FUNCTIONOBJECT_EXECUTE_HAS_NO_FORCE
    return functions().execute();
#else
#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
    if (forceWrite)
    {
        return functions().execute(); // there is no .write()-method
    }
    else
    {
        return functions().execute();
    }
#else
    return functions().execute(forceWrite);
#endif
#endif
}


bool functionObjectListProxy::start()
{
    Dbug<< this->name() << " functionObjectListProxy::start()" << endl;
    return functions().start();
}


bool functionObjectListProxy::end()
{
    Dbug<< this->name() << " functionObjectListProxy::end()" << endl;
    return functions().end();
}


bool functionObjectListProxy::read(const dictionary& dict)
{
    Dbug<< this->name() << " functionObjectListProxy::read()" << endl;
    return functions().read();
}


void functionObjectListProxy::writeSimple()
{}

} // namespace Foam

// ************************************************************************* //
