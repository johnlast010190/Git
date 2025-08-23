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
    2008-2011, 2013-2014, 2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "include/swak.H"

#include "loadTopoSetsFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "fvMesh/fvMesh.H"
#include "sets/topoSets/cellSet.H"
#include "sets/topoSets/faceSet.H"
#include "sets/topoSets/pointSet.H"
#include "include/swakTime.H"
#include "db/IOobjectList/IOobjectList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(loadTopoSetsFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        loadTopoSetsFunctionObject,
        dictionary
    );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

loadTopoSetsFunctionObject::loadTopoSetsFunctionObject
(
    const word &name,
    const Time& t,
    const dictionary& dict
)
:
    simpleFunctionObject(name,t,dict),
    forceLoading_(dict.lookup<bool>("forceLoading")),
    loadCellSets_(dict.lookup<bool>("loadCellSets")),
    loadFaceSets_(dict.lookup<bool>("loadFaceSets")),
    loadPointSets_(dict.lookup<bool>("loadPointSets")),
    writeSets_(dict.lookup<bool>("writeSets"))
{
}

bool loadTopoSetsFunctionObject::start()
{
    Info<< "Loading topoSets of " << obr().name()
        << " at start" << endl;

    loadSets();

    return true;
}

void loadTopoSetsFunctionObject::writeSimple()
{
    if (
        obr_.time().outputTime()
        &&
        writeSets_
    ) {
        writeAllSets<cellSet>(cellSetNames_);
        writeAllSets<faceSet>(faceSetNames_);
        writeAllSets<pointSet>(pointSetNames_);
    }
}

void loadTopoSetsFunctionObject::loadSets()
{
    if (loadCellSets_) {
        loadAllSets<cellSet>(cellSetNames_);
    }
    if (loadFaceSets_) {
        loadAllSets<faceSet>(faceSetNames_);
    }
    if (loadPointSets_) {
        loadAllSets<pointSet>(pointSetNames_);
    }
}

template<class TopoSetType>
void loadTopoSetsFunctionObject::writeAllSets(HashSet<word> &names)
{
    const polyMesh &mesh=dynamic_cast<const polyMesh&>(
        this->obr()
    );

    forAllConstIter(HashSet<word>,names,iter) {
        const word &name=iter.key();
        Info<< "Writing " << TopoSetType::typeName << " " << name << endl;

        TopoSetType &set=const_cast<TopoSetType &>(mesh.lookupObject<TopoSetType>(name));
        if (Pstream::parRun()) {
            set.instance()=
                (word("processor") + Foam::name(Pstream::myProcNo()))
                /
                mesh.time().timeName();
        } else {
            set.instance()=mesh.time().timeName();
        }
        set.write();
    }
}

template<class TopoSetType>
void loadTopoSetsFunctionObject::loadAllSets(HashSet<word> &names)
{
    Info<< "Loading sets of type " << TopoSetType::typeName << " from "
        << this->obr().name() << endl;

    const polyMesh &mesh=dynamic_cast<const polyMesh&>(
        this->obr()
    );

    // adapted from setsToZones.C

#ifdef FOAM_HAS_NO_FOUR_PARAMETER_FINDINSTANCE
    IOobjectList objects(mesh, mesh.facesInstance(), polyMesh::meshSubDir/"sets");
#else
    word setsInstance = mesh.time().findInstance
        (
            polyMesh::meshSubDir/"sets",
            word::null,
            IOobject::MUST_READ,
            mesh.facesInstance()
        );
    IOobjectList objects(mesh, setsInstance, polyMesh::meshSubDir/"sets");
#endif

    IOobjectList topoObjects(objects.lookupClass(TopoSetType::typeName));
    forAllConstIter(IOobjectList,topoObjects,iter) {
        const word &name=(*iter)->name();
        if (mesh.found(name)) {
            if (mesh.foundObject<TopoSetType>(name)) {
                if (forceLoading_) {
                    Info<< "Re-loading " << name << endl;
                    autoPtr<TopoSetType> set(new TopoSetType(*iter()));
                    Info<< " ... storing in " << mesh.name() <<
                        " (old version might still exist)" << endl;
                    set->store(set);
                    names.insert(name);
                } else {
                    WarningIn("loadTopoSetsFunctionObject::loadAllSets()")
                        << "There is already a " << TopoSetType::typeName
                            << " named " << name << " in " << mesh.name()
                            << ". Not loading"
                            << endl;
                }
            } else {
                FatalErrorIn("loadTopoSetsFunctionObject::loadAllSets()")
                    << "There is already an object named " << name
                        << " in " << mesh.name() << " but it is not a "
                        << TopoSetType::typeName
                        << endl
                        << exit(FatalError);
            }
        } else {
            Info<< "Reading " << name << endl;
            autoPtr<TopoSetType> set(new TopoSetType(*iter()));
            Info<< " ... storing in " << mesh.name() << endl;
            set->store(set);
            names.insert(name);
        }
    }
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

} // namespace Foam


// ************************************************************************* //
