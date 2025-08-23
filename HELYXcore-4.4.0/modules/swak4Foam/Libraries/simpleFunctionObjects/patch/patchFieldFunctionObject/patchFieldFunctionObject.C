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
    2008-2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "patchFieldFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchFieldFunctionObject, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

patchFieldFunctionObject::patchFieldFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    patchFunctionObject(name,t,dict),
    fieldNames_(0),
    scalarFields_(0),
    vectorFields_(0),
    sphericalTensorFields_(0),
    symmTensorFields_(0),
    tensorFields_(0)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void patchFieldFunctionObject::writeSimple()
{
    processAndWrite<scalar>(scalarFields_);
    processAndWrite<vector>(vectorFields_);
    processAndWrite<sphericalTensor>(sphericalTensorFields_);
    processAndWrite<symmTensor>(symmTensorFields_);
    processAndWrite<tensor>(tensorFields_);
}

bool patchFieldFunctionObject::start()
{
    fieldNames_ = wordList(dict_.lookup("fields"));

    boolList foundFields(fieldNames_.size(), false);
    findFields<volScalarField>(scalarFields_, foundFields);
    findFields<volVectorField>(vectorFields_, foundFields);
    findFields<volSphericalTensorField>(sphericalTensorFields_, foundFields);
    findFields<volSymmTensorField>(symmTensorFields_, foundFields);
    findFields<volTensorField>(tensorFields_, foundFields);

    label validFieldi=0;
    forAll(fieldNames_, fieldi)
    {
        if (foundFields[fieldi])
        {
            fieldNames_[validFieldi++] = fieldNames_[fieldi];
        }
        else
        {
            WarningIn("probes::read()")
                << "Unknown field " << fieldNames_[fieldi]
                << " when reading dictionary " << dict_.name() << endl
                << "    Can only probe registered volScalar, volVector, "
                << "volSphericalTensor, volSymmTensor and volTensor fields"
                << nl << endl;
        }
    }

    fieldNames_.setSize(validFieldi);


    if (Pstream::master())
    {
        if (debug)
        {
            Pout<< "Probing fields:" << fieldNames_ << nl;
        }

        // Check if any fieldNames have been removed. If so close
        // the file.
    }

    patchFunctionObject::start();

    return true;
}

wordList patchFieldFunctionObject::fileNames()
{
    return fieldNames_;
}

} // namespace Foam

// ************************************************************************* //
