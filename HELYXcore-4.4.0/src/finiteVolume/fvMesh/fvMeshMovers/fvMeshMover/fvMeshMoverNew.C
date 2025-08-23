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
    (c) 2021-2022 OpenFOAM Foundation
    (c) 2023-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMesh/fvMeshMovers/none/fvMeshMoversNone.H"
#include "db/dynamicLibrary/dlLibraryTable/dlLibraryTable.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fvMeshMover> Foam::fvMeshMover::New(fvMesh& mesh)
{
    IOobject dictHeader
    (
        "dynamicMeshDict",
        mesh.time().constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    if (dictHeader.typeHeaderOk<IOdictionary>())
    {
        IOdictionary dict(dictHeader);

        // Check if the old keyword 'dynamicFvMesh' is present
        if (dict.found("dynamicFvMesh"))
        {
            FatalErrorInFunction
                << "Old keyword 'dynamicFvMesh' found in the dynamicMeshDict "
                << "dictionary. Please, make sure that the new dynamic mesh "
                << "structure is used and remove that keyword!"
                << exit(FatalError);
        }

        if (dict.found("mover"))
        {
            const dictionary& moverDict = dict.subDict("mover");

            const word fvMeshMoverTypeName(moverDict.lookup("type"));

            Info<< "Selecting fvMeshMover " << fvMeshMoverTypeName << endl;

            libs.open(moverDict, "libs");

            auto ctor =
                ctorTableLookup
                (
                    "fvMeshMover",
                    fvMeshConstructorTable_(),
                    fvMeshMoverTypeName
                );

            return autoPtr<fvMeshMover>(ctor(mesh));
        }
    }

    return autoPtr<fvMeshMover>(new fvMeshMovers::none(mesh));
}


// ************************************************************************* //
