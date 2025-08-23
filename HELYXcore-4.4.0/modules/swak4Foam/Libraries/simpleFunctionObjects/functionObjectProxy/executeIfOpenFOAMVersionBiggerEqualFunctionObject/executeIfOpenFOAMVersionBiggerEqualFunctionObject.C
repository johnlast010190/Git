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
    2011, 2013, 2015-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "executeIfOpenFOAMVersionBiggerEqualFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "meshes/polyMesh/polyMesh.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "include/swakTime.H"

#include "foamVersion4swak.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(executeIfOpenFOAMVersionBiggerEqualFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        executeIfOpenFOAMVersionBiggerEqualFunctionObject,
        dictionary
    );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

executeIfOpenFOAMVersionBiggerEqualFunctionObject::executeIfOpenFOAMVersionBiggerEqualFunctionObject
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
    // do it here to avoid the superclass-read being read twice
    readData(dict);

#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
    start();
#endif
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label executeIfOpenFOAMVersionBiggerEqualFunctionObject::toLabel(const string &v) {
    if (v=="") {
        return -1;
    } else if (v=="x") {
        return labelMax;
    } else {
        IStringStream i(v);

        return readLabel(i);
    }
}

bool executeIfOpenFOAMVersionBiggerEqualFunctionObject::condition()
{
#define TOSTRING(x) string(#x)

    label foamVersionPatch=-1;
    if (isdigit(TOSTRING(FOAM_VERSION4SWAK_PATCH)[0])) {
        foamVersionPatch=toLabel(TOSTRING(FOAM_VERSION4SWAK_PATCH));
    }
    if (majorVersion_>FOAM_VERSION4SWAK_MAJOR) {
        return false;
    } else if (majorVersion_<FOAM_VERSION4SWAK_MAJOR) {
        return true;
    } else if (minorVersion_>FOAM_VERSION4SWAK_MINOR) {
        return false;
    } else if (minorVersion_<FOAM_VERSION4SWAK_MINOR) {
        return true;
    } else if (foamVersionPatch<0) {
        return true;
    } else {
        return patchVersion_<=foamVersionPatch;
    }
#undef TOSTRING
}

void executeIfOpenFOAMVersionBiggerEqualFunctionObject::readData(const dictionary& dict)
{
    majorVersion_=toLabel(string(dict.lookup("majorVersion")));
    minorVersion_=toLabel(string(dict.lookup("minorVersion")));
    patchVersion_=toLabel(dict.lookupOrDefault<string>("patchVersion",""));
}

bool executeIfOpenFOAMVersionBiggerEqualFunctionObject::read(const dictionary& dict)
{
    readData(dict);
    return conditionalFunctionObjectListProxy::read(dict);
}

} // namespace Foam

// ************************************************************************* //
