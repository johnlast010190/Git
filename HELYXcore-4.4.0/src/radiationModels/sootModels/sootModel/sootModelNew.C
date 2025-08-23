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
    (c) 2013-2020 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "db/error/error.H"
#include "sootModels/sootModel/sootModel.H"
#include "sootModels/noSoot/noSoot.H"
#include "basicThermo/basicThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::radiationModels::sootModel>
Foam::radiationModels::sootModel::New
(
    const dictionary& dict,
    const fvMesh& mesh
)
{
    // Get the soot model type name
    const word modelType =
        dict.lookupOrDefault<word>
        (
            sootModel::typeName,
            sootModels::noSoot::typeName
        );
    Info<< "Selecting soot model " << modelType << endl;

    // Lookup both possible model names
    const auto ctor =
        ctorTableLookup
        (
            sootModel::typeName + " type",
            dictionaryConstructorTable_(),
            modelType
        );
    return autoPtr<sootModel>(ctor(dict, mesh, modelType));
}


// ************************************************************************* //
