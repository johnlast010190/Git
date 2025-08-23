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
    (c) 2017-2018 OpenFOAM Foundation
    (c) 2019 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "global/fileOperations/fileOperationInitialise/fileOperationInitialise.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(fileOperationInitialise, 0);
    defineRunTimeSelectionTable(fileOperationInitialise, word);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::fileOperationInitialise::fileOperationInitialise
(
    int& argc,
    char**& argv
)
{}


Foam::autoPtr<Foam::fileOperations::fileOperationInitialise>
Foam::fileOperations::fileOperationInitialise::New
(
    const word& type,
    int& argc,
    char**& argv
)
{
    DebugInFunction << "Constructing fileOperationInitialise" << endl;

    const auto ctor =
        ctorTableLookup
        (
            "fileOperationInitialise type",
            wordConstructorTable_(),
            type
        );
    return autoPtr<fileOperationInitialise>(ctor(argc, argv));
}


// ************************************************************************* //
