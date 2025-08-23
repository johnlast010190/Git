/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2015 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/sensor/patchInternalFieldValue/patchInternalFieldValue.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sensorTypes
{
    defineTypeNameAndDebug(patchInternalFieldValue, 0);

    sensor<scalar>::addmeshConstructorToTable<sensorTypes::patchInternalFieldValue>
        addpatchInternalFieldValueMeshConstructorToTable_;
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sensorTypes::patchInternalFieldValue::patchInternalFieldValue
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    sensor<scalar>(mesh, dict),
    patchName_()
{
    Istream& is(dict.lookup("patchName"));
    is  >> patchName_;
}


Foam::sensorTypes::patchInternalFieldValue::patchInternalFieldValue(const patchInternalFieldValue& cv)
:
    sensor<scalar>(cv),
    patchName_(cv.patchName_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sensorTypes::patchInternalFieldValue::~patchInternalFieldValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::sensorTypes::patchInternalFieldValue::valueField() const
{
    const volScalarField& field(mesh_.lookupObject<volScalarField>(fieldName()));
    label patchID = mesh_.boundaryMesh().findPatchID(patchName_);

    return field.boundaryField()[patchID].patchInternalField();
}


// ************************************************************************* //
