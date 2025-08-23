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
    (c) 2016 OpenCFD Ltd.
    (c) 2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "helpSolver.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace helpTypes
    {
        defineTypeNameAndDebug(helpSolver, 0);
        addNamedToRunTimeSelectionTable
        (
            helpType,
            helpSolver,
            dictionary,
            solver
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::helpTypes::helpSolver::helpSolver()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::helpTypes::helpSolver::~helpSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::helpTypes::helpSolver::init()
{
    helpType::init();

    argList::addBoolOption
    (
        "read",
        "read solver type from the system/controlDict"
    );
}


void Foam::helpTypes::helpSolver::execute
(
    const argList& args,
    const fvMesh& mesh
)
{
    word solver(word::null);

    if (args.optionReadIfPresent("browse", solver))
    {
        displayDoc(solver, ".*solvers/.*Foam/", true, "C");
    }
    else if (args.optionFound("read"))
    {
        solver = mesh.time().controlDict().lookup<word>("application");
        displayDoc(solver, ".*solvers/.*Foam/", true, "C");
    }
    else
    {
        displayDocOptions(".*solvers/.*Foam/", true, "C");
    }
}


// ************************************************************************* //
