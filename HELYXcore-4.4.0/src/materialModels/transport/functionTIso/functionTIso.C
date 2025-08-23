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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2019-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "functionTIso.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(functionTIso, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        functionTIso,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionTIso::functionTIso
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
:
    materialModel(obr, dict, phaseName, specieName, name)
{
    functionTIso::read();
}


Foam::autoPtr<Foam::functionTIso>
Foam::functionTIso::clone() const
{
    return autoPtr<functionTIso>(new functionTIso(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::baseModels<Foam::scalar>* Foam::functionTIso::castScalarModel
(
    const word& modelName
)
{
    castMaterial(modelName, kappa)

    return nullptr;
}


Foam::baseModels<Foam::vector>* Foam::functionTIso::castVectorModel
(
    const word& modelName
)
{
    castMaterial(modelName, vKappa)

    return nullptr;
}


bool Foam::functionTIso::read()
{
    initialise("T");
    kappa_.reset(Function1<scalar>::New("kappa", *dict_));

    return true;
}


// ************************************************************************* //
