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
    2011, 2013, 2015-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "executeIfObjectExistsFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "meshes/polyMesh/polyMesh.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "include/swakTime.H"

#include "db/objectRegistry/objectRegistry.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(executeIfObjectExistsFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        executeIfObjectExistsFunctionObject,
        dictionary
    );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

executeIfObjectExistsFunctionObject::executeIfObjectExistsFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    conditionalFunctionObjectListProxy(
        name,
        t,
        dict
    )
{
    readParameters(dict);

#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
    start();
#endif
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool executeIfObjectExistsFunctionObject::condition()
{
    if (writeDebug()) {
        Info<< "Looking for object " << objectName_
            << " -> " << obr().foundObject<IOobject>(objectName_) << endl;
    }
    if (! obr().foundObject<IOobject>(objectName_)) {
        return ! objectShouldExist_;
    } else if (checkType_) {
        const IOobject &theOb=obr().lookupObject<IOobject>(objectName_);
        if (writeDebug()) {
            Info<< "Type of " << objectName_ << " is "
                << theOb.type()
                << ". Looking for " << objectType_ << endl;
        }

        if (theOb.type()==objectType_) {
            return objectShouldExist_;
        } else {
            return ! objectShouldExist_;
        }
    } else {
        return objectShouldExist_;
    }
}

bool executeIfObjectExistsFunctionObject::read(const dictionary& dict)
{
    readParameters(dict);
    return conditionalFunctionObjectListProxy::read(dict);
}

void executeIfObjectExistsFunctionObject::readParameters(const dictionary &dict)
{
    objectName_=word(dict.lookup("objectName"));
    checkType_=dict.lookup<bool>("checkType");
    if (checkType_) {
        objectType_=word(dict.lookup("objectType"));
    } else {
        objectType_="";
    }
    objectShouldExist_=dict.lookup<bool>("objectShouldExist");
}

} // namespace Foam

// ************************************************************************* //
