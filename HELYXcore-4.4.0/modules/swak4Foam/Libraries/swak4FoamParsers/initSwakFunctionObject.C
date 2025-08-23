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
    2008-2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "initSwakFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "meshes/polyMesh/polyMesh.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "include/swakTime.H"

#include "CommonValueExpressionDriver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(initSwakFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        initSwakFunctionObject,
        dictionary
    );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

initSwakFunctionObject::initSwakFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    dict_(dict)
{
    word regionName=
        dict.lookupOrDefault<word>("region",polyMesh::defaultRegion);

    const fvMesh &mesh=dynamic_cast<const fvMesh &>(
        t.lookupObject<objectRegistry>(regionName)
    );

    CommonValueExpressionDriver::resetDefaultMesh(mesh);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool initSwakFunctionObject::start()
{
    return true;
}

bool initSwakFunctionObject::execute(const bool forceWrite)
{
    return true;
}

#if defined(FOAM_FUNCTIONOBJECT_EXECUTE_HAS_NO_FORCE) || defined(FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START)
bool initSwakFunctionObject::execute()
{
    return execute(false);
}
#endif

#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
bool initSwakFunctionObject::write()
{
    return execute(true);
}
#endif

bool initSwakFunctionObject::read(const dictionary& dict)
{
    return true;
}

} // namespace Foam

// ************************************************************************* //
