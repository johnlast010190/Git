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
    (c) 2021-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "kappaSutherland.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(kappaSutherland, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        kappaSutherland,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kappaSutherland::kappaSutherland
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
    kappaSutherland::read();
}


Foam::autoPtr<Foam::kappaSutherland>
Foam::kappaSutherland::clone() const
{
    return autoPtr<kappaSutherland>(new kappaSutherland(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::kappaSutherland::updateTable(const word& modelName)
{
    // Required dependencies
    const wordList depNames
    ({
        CvModel::typeName,
        muModel::typeName,
        RModel::typeName
    });

    // Dependency indices
    const labelList depInds({Cv, mu, R});

    // Create links to dependent models
    fill(modelName, kappaModel::typeName, depNames, depInds, 0);
}


Foam::baseModels<Foam::scalar>* Foam::kappaSutherland::castScalarModel
(
    const word& modelName
)
{
    castMaterial(modelName, kappa)

    return nullptr;
}


bool Foam::kappaSutherland::read()
{
    return true;
}


// ************************************************************************* //
