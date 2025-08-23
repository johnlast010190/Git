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
    2011-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "executeIfSetOrZoneFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "meshes/polyMesh/zones/cellZone/cellZone.H"
#include "meshes/polyMesh/zones/faceZone/faceZone.H"
#include "meshes/polyMesh/zones/pointZone/pointZone.H"
#include "sets/topoSets/cellSet.H"
#include "sets/topoSets/faceSet.H"
#include "sets/topoSets/pointSet.H"
#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(executeIfSetOrZoneFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        executeIfSetOrZoneFunctionObject,
        dictionary
    );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

executeIfSetOrZoneFunctionObject::executeIfSetOrZoneFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    conditionalFunctionObjectListProxy(name, t, dict),
    mesh_(dynamicCast<const polyMesh&>(obr())),
    loadAndCacheMissingSets_(dict.lookup<bool>("loadAndCacheMissingSets"))
{
    // do it here to avoid the superclass-read being read twice
    readSetsAndZones(dict);

#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
    start();
#endif
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
bool executeIfSetOrZoneFunctionObject::hasSet(const word &name)
{
    if (mesh_.foundObject<T>(name))
    {
        Dbug<< name << " of type " << T::typeName
            << " already in memory" << endl;
        return true;
    }
    else
    {
        Dbug<< name << " of type " << T::typeName << " not in memory" << endl;
        if (loadAndCacheMissingSets_)
        {
            Dbug<< "Loading " << name << endl;
            autoPtr<T> s(new T(mesh_, name, IOobject::READ_IF_PRESENT));

            if (s->headerOk())
            {
                Dbug<< "Storing " << name << " in mesh" << endl;
                s->store(s);
                return true;
            }
            else
            {
                Dbug<< "No valid " << name << endl;
            }
        }
        return false;
    }
}


bool executeIfSetOrZoneFunctionObject::condition()
{
    forAllConstIter(dictionary, setsAndZones_, it)
    {
        const entry& e = *it;
        const word& name = e.keyword();
        word typ(e.stream());

        Dbug<< "Typ: " << typ << " Name: " << name << endl;

    if (typ == "cellZone")
    {
#ifdef FOAM_ZONEMESH_HAS_NO_FINDINDEX
        if (mesh_.cellZones().findZoneID(name) < 0)
        {
#else
        if (mesh_.cellZones().findIndex(name) < 0)
        {
#endif
            Dbug<< "No " << name << " of type " << typ << endl;
            return false;
        }
    }
    else if (typ == "faceZone")
    {
#ifdef FOAM_ZONEMESH_HAS_NO_FINDINDEX
        if (mesh_.faceZones().findZoneID(name) < 0)
        {
#else
        if (mesh_.faceZones().findIndex(name) < 0)
        {
#endif
            Dbug<< "No " << name << " of type " << typ << endl;
            return false;
        }
    }
    else if (typ == "pointZone")
    {
#ifdef FOAM_ZONEMESH_HAS_NO_FINDINDEX
        if (mesh_.pointZones().findZoneID(name) < 0)
        {
#else
        if (mesh_.pointZones().findIndex(name) < 0)
        {
#endif
                Dbug<< "No " << name << " of type " << typ << endl;
                return false;
            }
        }
        else if (typ == "cellSet")
        {
            if (!this->hasSet<cellSet>(name))
            {
                Dbug<< "No " << name << " of type " << typ << endl;
                return false;
            }
        }
        else if (typ == "faceSet")
        {
            if (!this->hasSet<faceSet>(name))
            {
                Dbug<< "No " << name << " of type " << typ << endl;
                return false;
            }
        }
        else if (typ == "pointSet")
        {
            if (!this->hasSet<pointSet>(name))
            {
                Dbug<< "No " << name << " of type " << typ << endl;
                return false;
            }
        }
        else
        {
            FatalErrorIn("executeIfSetOrZoneFunctionObject::condition()")
                << "Unimplemented type " << typ << " for " << name << endl
                << "Implemented are faceZone, faceSet, cellZone,"
                << " cellSet, pointZone, pointSet" << endl
                << abort(FatalError);
        }
    }
    return true;
}


void executeIfSetOrZoneFunctionObject::readSetsAndZones(const dictionary& dict)
{
    setsAndZones_ = dict.subDict("setsAndZones");
    loadAndCacheMissingSets_=
        dict.lookup<bool>("loadAndCacheMissingSets");
}


bool executeIfSetOrZoneFunctionObject::read(const dictionary& dict)
{
    readSetsAndZones(dict);
    return conditionalFunctionObjectListProxy::read(dict);
}

} // namespace Foam

// ************************************************************************* //
