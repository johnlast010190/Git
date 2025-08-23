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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "matrices/lduMatrix/lduMatrix/lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineRunTimeSelectionTable(lduMatrix::preconditioner, symMatrix);
    defineRunTimeSelectionTable(lduMatrix::preconditioner, asymMatrix);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word Foam::lduMatrix::preconditioner::getName
(
    const dictionary& solverControls
)
{
    word name;

    // handle primitive or dictionary entry
    const entry& e = solverControls.lookupEntry("preconditioner", false, false);
    if (e.isDict())
    {
        name = e.dict().lookup<word>("preconditioner");
    }
    else
    {
        e.stream() >> name;
    }

    return name;
}


Foam::autoPtr<Foam::lduMatrix::preconditioner>
Foam::lduMatrix::preconditioner::New
(
    const solver& sol,
    const dictionary& solverControls
)
{
    word name;

    // handle primitive or dictionary entry
    const entry& e = solverControls.lookupEntry("preconditioner", false, false);
    if (e.isDict())
    {
        name = e.dict().lookup<word>("preconditioner");
    }
    else
    {
        e.stream() >> name;
    }

    const dictionary& controls = e.isDict() ? e.dict() : dictionary::null;

    if (sol.matrix().symmetric())
    {
        auto ctor =
            ctorTableLookup
            (
                "symmetric matrix preconditioner",
                symMatrixConstructorTable_(),
                name,
                controls
            );

        return autoPtr<lduMatrix::preconditioner>
        (
            ctor
            (
                sol,
                controls
            )
        );
    }
    else if (sol.matrix().asymmetric())
    {
        auto ctor =
            ctorTableLookup
            (
                "asymmetric matrix preconditioner",
                asymMatrixConstructorTable_(),
                name,
                controls
            );

        return autoPtr<lduMatrix::preconditioner>
        (
            ctor
            (
                sol,
                controls
            )
        );
    }
    else
    {
        FatalIOErrorInFunction
        (
            controls
        )   << "cannot solve incomplete matrix, "
               "no diagonal or off-diagonal coefficient"
            << exit(FatalIOError);

        return autoPtr<lduMatrix::preconditioner>(nullptr);
    }
}


// ************************************************************************* //
