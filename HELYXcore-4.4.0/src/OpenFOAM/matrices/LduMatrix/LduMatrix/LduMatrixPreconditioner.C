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

\*---------------------------------------------------------------------------*/

#include "matrices/LduMatrix/LduMatrix/LduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::autoPtr<typename Foam::LduMatrix<Type, DType, LUType>::preconditioner>
Foam::LduMatrix<Type, DType, LUType>::preconditioner::New
(
    const solver& sol,
    const dictionary& preconditionerDict
)
{
    word preconditionerName = preconditionerDict.lookup("preconditioner");

    if (sol.matrix().symmetric())
    {
        typename symMatrixConstructorTable::iterator constructorIter =
            symMatrixConstructorTable_().find(preconditionerName);

        if (constructorIter == symMatrixConstructorTable_().end())
        {
            FatalIOErrorInFunction
            (
                preconditionerDict
            )   << "Unknown symmetric matrix preconditioner "
                << preconditionerName << endl << endl
                << exit(FatalIOError);
        }

        return autoPtr<typename LduMatrix<Type, DType, LUType>::preconditioner>
        (
            constructorIter->second
            (
                sol,
                preconditionerDict
            )
        );
    }
    else if (sol.matrix().asymmetric())
    {
        typename asymMatrixConstructorTable::iterator constructorIter =
            asymMatrixConstructorTable_().find(preconditionerName);

        if (constructorIter == asymMatrixConstructorTable_().end())
        {
            FatalIOErrorInFunction
            (
                preconditionerDict
            )   << "Unknown asymmetric matrix preconditioner "
                << preconditionerName << endl << endl
                << exit(FatalIOError);
        }

        return autoPtr<typename LduMatrix<Type, DType, LUType>::preconditioner>
        (
            constructorIter->second
            (
                sol,
                preconditionerDict
            )
        );
    }
    else
    {
        FatalIOErrorInFunction
        (
            preconditionerDict
        )   << "cannot preconditione incomplete matrix, "
               "no diagonal or off-diagonal coefficient"
            << exit(FatalIOError);

        return autoPtr<typename LduMatrix<Type, DType, LUType>::preconditioner>
        (
            nullptr
        );
    }
}


// ************************************************************************* //
