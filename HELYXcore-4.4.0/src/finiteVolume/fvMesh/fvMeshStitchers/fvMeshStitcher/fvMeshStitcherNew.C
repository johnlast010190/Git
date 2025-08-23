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
    (c) 2022 OpenFOAM Foundation
    (c) 2022-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMesh/fvMeshStitchers/fvMeshStitcher/fvMeshStitcher.H"
#include "db/dynamicLibrary/dlLibraryTable/dlLibraryTable.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fvMeshStitcher> Foam::fvMeshStitcher::New
(
    fvMesh& mesh,
    const bool changers
)
{
    // Determine if the mesh is actually changing and load the
    // fvMeshStitchers library if so
    bool changing = changers;
    if (changers)
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

        changing = changers && dictHeader.typeHeaderOk<IOdictionary>();

        if (changing)
        {
            libs.open("lib" + typeName + "s.so");
        }
    }

    // Select a non-changing or changing stitcher as appropriate
    for (const auto& p : fvMeshConstructorTable_())
    {
        autoPtr<fvMeshStitcher> stitcherPtr(p.second(mesh));

        if (stitcherPtr->changing() == changing)
        {
            return stitcherPtr;
        }
    }

    // Error if an appropriate stitcher was not found
    FatalErrorInFunction
        << typeName << " for " << (changing ? "" : "non-")
        << "changing mesh not found " << nl << nl
        << exit(FatalError);

    return autoPtr<fvMeshStitcher>(nullptr);
}


// ************************************************************************* //
