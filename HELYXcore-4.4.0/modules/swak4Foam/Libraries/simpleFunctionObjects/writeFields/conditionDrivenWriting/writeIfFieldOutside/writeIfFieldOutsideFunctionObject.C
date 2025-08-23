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

 SWAK Revision: $Id: writeIfFieldOutsideFunctionObject.C,v 78b0d113b99b 2013-02-25 16:12:41Z bgschaid $
\*---------------------------------------------------------------------------*/

#include "writeIfFieldOutsideFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "meshes/polyMesh/polyMesh.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "include/swakTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(writeIfFieldOutsideFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        writeIfFieldOutsideFunctionObject,
        dictionary
    );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

writeIfFieldOutsideFunctionObject::writeIfFieldOutsideFunctionObject
(
    const word &name,
    const Time& t,
    const dictionary& dict
)
:
    conditionDrivenWritingFunctionObject(name,t,dict),
    fieldName_(""),
    maximum_(pTraits<scalar>::max),
    minimum_(pTraits<scalar>::min)
{
}

bool writeIfFieldOutsideFunctionObject::start()
{
    conditionDrivenWritingFunctionObject::start();

    fieldName_=word(dict_.lookup("fieldName"));
    minimum_=dict_.lookup<scalar>("minimum");
    maximum_=dict_.lookup<scalar>("maximum");

    Info<< "Checking for field " << fieldName_ << " in range [ " << minimum_
        << " , " << maximum_ << " ] " << endl;

    return true;
}

bool writeIfFieldOutsideFunctionObject::checkStartWriting()
{
    bool outside=
        check<volScalarField>() ||
        check<volVectorField>() ||
        check<volSphericalTensorField>() ||
        check<volSymmTensorField>() ||
        check<volTensorField>();

    reduce(outside,orOp<bool>());

    return outside;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

} // namespace Foam

// ************************************************************************* //
