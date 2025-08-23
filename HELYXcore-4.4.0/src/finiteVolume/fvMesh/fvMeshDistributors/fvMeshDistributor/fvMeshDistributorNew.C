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
    (c) 2021-2023 OpenFOAM Foundation
    (c) 2023-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMesh/fvMeshDistributors/none/fvMeshDistributorsNone.H"
#include "db/dynamicLibrary/dlLibraryTable/dlLibraryTable.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fvMeshDistributor> Foam::fvMeshDistributor::New
(
    fvMesh& mesh
)
{
    // Only construct a real distributor when running in parallel
    // otherwise return a fvMeshDistributors::none
    if (Pstream::parRun())
    {
        //typeIOobject<IOdictionary> dictHeader
        //(
            IOobject dictHeader
            (
                "dynamicMeshDict",
                mesh.time().constant(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            );
        //);

        if (dictHeader.headerOk())
        {
            IOdictionary dict(dictHeader);

            if (dict.found("distributor"))
            {
                const dictionary& distributorDict = dict.subDict("distributor");

                const word fvMeshDistributorTypeName
                (
                    distributorDict.lookup("type")
                );

                Info<< "Selecting fvMeshDistributor "
                    << fvMeshDistributorTypeName << endl;

                libs.open(distributorDict, "libs");

                auto ctor = ctorTableLookup("fvMeshDistributor", fvMeshConstructorTable_(), fvMeshDistributorTypeName);
                return autoPtr<fvMeshDistributor>(ctor(mesh));
            }
        }
    }

    return autoPtr<fvMeshDistributor>(new fvMeshDistributors::none(mesh));
}


// ************************************************************************* //
