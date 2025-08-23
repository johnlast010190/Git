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

Contributors/Copyright:
    2011-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "executeIfPatchFunctionObject.H"
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
    defineTypeNameAndDebug(executeIfPatchFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        executeIfPatchFunctionObject,
        dictionary
    );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

executeIfPatchFunctionObject::executeIfPatchFunctionObject
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
    ),
    mesh_(
    dynamicCast<const polyMesh&>(obr())
    )
{
    // do it here to avoid the superclass-read being read twice
    readPatches(dict);

#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
    start();
#endif
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool executeIfPatchFunctionObject::condition()
{
    forAll(patchNames_,i) {
        if (mesh_.boundaryMesh().findPatchID(patchNames_[i])<0) {
            return false;
        }
    }
    return true;
}

void executeIfPatchFunctionObject::readPatches(const dictionary& dict)
{
    patchNames_=wordList(dict.lookup("patchNames"));
}

bool executeIfPatchFunctionObject::read(const dictionary& dict)
{
    readPatches(dict);
    return conditionalFunctionObjectListProxy::read(dict);
}

} // namespace Foam

// ************************************************************************* //
