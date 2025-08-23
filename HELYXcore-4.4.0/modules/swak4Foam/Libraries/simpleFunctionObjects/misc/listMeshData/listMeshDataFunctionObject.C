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
    2008-2011, 2013-2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "include/swak.H"
#ifdef FOAM_FVMESH_IS_DATAOBJECT

#include "listMeshDataFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "meshes/polyMesh/polyMesh.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "include/swakTime.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(listMeshDataFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        listMeshDataFunctionObject,
        dictionary
    );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

listMeshDataFunctionObject::listMeshDataFunctionObject
(
    const word &name,
    const Time& t,
    const dictionary& dict
)
:
    simpleFunctionObject(name,t,dict)
{
}

bool listMeshDataFunctionObject::start()
{
    Info<< "Content of mesh data " << obr().name()
        << " at start" << endl;

    dumpData();

    return true;
}

void listMeshDataFunctionObject::writeSimple()
{
    Info<< "Content of mesh data " << obr().name()
        << endl;

    dumpData();
}

void listMeshDataFunctionObject::dumpData()
{
    const data &theData=dynamicCast<const fvMesh&>(obr());

    Info<< theData << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

} // namespace Foam

#endif

// ************************************************************************* //
