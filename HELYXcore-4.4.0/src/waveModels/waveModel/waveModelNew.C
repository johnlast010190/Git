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
    (c) 2015 IH-Cantabria
    (c) 2016 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "waveModel/waveModel.H"
#include "fvMesh/fvMesh.H"

Foam::autoPtr<Foam::waveModel> Foam::waveModel::New
(
    const word& dictName,
    const fvMesh& mesh,
    const polyPatch& patch
)
{
    IOdictionary waveDict
    (
        IOobject
        (
            dictName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false                   // Not registering
        )
    );

    word modelType = "none";
    dictionary patchDict;
    if (waveDict.found(patch.name()))
    {
        patchDict = waveDict.subDict(patch.name());
        modelType = patchDict.lookup<word>("waveModel");
    }
    else
    {
        FatalIOErrorInFunction(waveDict)
            << "Dictionary entry for patch " << patch.name() << " not found"
            << exit(FatalIOError);
    }

    Info<< "Selecting waveModel " << modelType << endl;

    const auto ctor =
        ctorTableLookup
        (
            "waveModel type",
            patchConstructorTable_(),
            modelType
        );
    return autoPtr<waveModel>(ctor(patchDict, mesh, patch));
}


Foam::tmp<Foam::waveModel> Foam::waveModel::lookupOrCreate
(
    const polyPatch& patch,
    const fvMesh& mesh,
    const word& waveDictName
)
{
    const word modelName = waveModel::modelName(patch.name());

    if (!mesh.foundObject<waveModel>(modelName))
    {
        autoPtr<waveModel> model(waveModel::New(waveDictName, mesh, patch));
        waveModel* waveModelPtr = model.ptr();
        waveModelPtr->store();
        waveModelPtr->info(Info);
    }

    return mesh.lookupObject<waveModel>(modelName);
}


// ************************************************************************* //
