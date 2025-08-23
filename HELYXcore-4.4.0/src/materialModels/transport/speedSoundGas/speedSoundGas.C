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
    (c) 2021-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "speedSoundGas.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(speedSoundGas, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        speedSoundGas,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::speedSoundGas::speedSoundGas
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
    sMod_.setSize(modelsEnumSize_);
    dep_.setSize(1);
    speedSoundGas::read();
}


Foam::autoPtr<Foam::speedSoundGas>
Foam::speedSoundGas::clone() const
{
    return autoPtr<speedSoundGas>(new speedSoundGas(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::speedSoundGas::updateTable(const word& modelName)
{
    // Required dependencies
    const wordList depNames
    ({
        gammaModel::typeName,
        RModel::typeName
    });

    // Dependency indices
    const labelList depInds({gamma, R});

    // Create links to dependent models
    fill(modelName, c0Model::typeName, depNames, depInds, 0);
}


Foam::baseModels<Foam::scalar>* Foam::speedSoundGas::castScalarModel
(
    const word& modelName
)
{
    castMaterial(modelName, c0)

    return nullptr;
}


bool Foam::speedSoundGas::read()
{
    initialise("T");

    return true;
}


// ************************************************************************* //
