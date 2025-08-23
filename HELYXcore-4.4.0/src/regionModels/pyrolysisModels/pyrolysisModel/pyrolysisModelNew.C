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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "pyrolysisModels/pyrolysisModel/pyrolysisModel.H"
#include "fvMesh/fvMesh.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::regionModels::pyrolysisModel>
Foam::regionModels::pyrolysisModel::New
(
    const fvMesh& mesh,
    const word& regionType
)
{
    // get model name, but do not register the dictionary
    const word modelType
    (
        IOdictionary
        (
            IOobject
            (
                regionType + "Properties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("pyrolysisModel")
    );
    if (modelType == "none")
    {
        return autoPtr<pyrolysisModel>
        (
            new pyrolysisModel(modelType, mesh, regionType)
        );
    }

    Info<< "Selecting pyrolysisModel " << modelType << endl;

    const auto ctor =
        ctorTableLookup
        (
            "pyrolysisModel type",
            meshConstructorTable_(),
            modelType
        );
    return autoPtr<pyrolysisModel>(ctor(modelType, mesh, regionType));
}


Foam::autoPtr<Foam::regionModels::pyrolysisModel>
Foam::regionModels::pyrolysisModel::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& regionType
)
{
    const word modelType = dict.lookup("pyrolysisModel");
    if (modelType == "none")
    {
        return autoPtr<pyrolysisModel>
        (
            new pyrolysisModel(modelType, mesh, dict, regionType)
        );
    }

    Info<< "Selecting pyrolysisModel " << modelType << endl;

    const auto ctor =
        ctorTableLookup
        (
            "pyrolysisModel type",
            dictionaryConstructorTable_(),
            modelType
        );
    return autoPtr<pyrolysisModel>(ctor(modelType, mesh, dict, regionType));
}


// ************************************************************************* //
