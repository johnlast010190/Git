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
    2012-2013, 2015-2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "writeAndEndFieldRangeFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "meshes/polyMesh/polyMesh.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "include/swakTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(writeAndEndFieldRangeFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        writeAndEndFieldRangeFunctionObject,
        dictionary
    );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

writeAndEndFieldRangeFunctionObject::writeAndEndFieldRangeFunctionObject
(
    const word &name,
    const Time& t,
    const dictionary& dict
)
:
    writeAndEndFunctionObject(name,t,dict),
    fieldName_(""),
    maximum_(pTraits<scalar>::max),
    minimum_(pTraits<scalar>::min)
{
}

bool writeAndEndFieldRangeFunctionObject::start()
{
    writeAndEndFunctionObject::start();

    fieldName_=word(dict_.lookup("fieldName"));
    minimum_=dict_.lookup<scalar>("minimum");
    maximum_=dict_.lookup<scalar>("maximum");

    Info<< "Checking for field " << fieldName_ << " in range [ " << minimum_
        << " , " << maximum_ << " ] " << endl;

    return true;
}

bool writeAndEndFieldRangeFunctionObject::endRunNow()
{
    if (
        check<volScalarField>()
        ||
        check<volVectorField>()
        ||
        check<volSphericalTensorField>()
        ||
        check<volSymmTensorField>()
        ||
        check<volTensorField>()
    ) {
        return true;
    } else {
        return false;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

} // namespace Foam

// ************************************************************************* //
