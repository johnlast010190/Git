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
    (c) 2004-6 H. Jasak All rights reserved
    (c) 2024 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "matrices/blockLduMatrix/BlockLduPrecons/BlockLduPrecon/BlockLduPrecon.H"
#include "matrices/blockLduMatrix/BlockLduPrecons/BlockNoPrecon/blockNoPrecons.H"

template<class Type>
class BlockNoPrecon;

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::BlockLduPrecon<Type>> Foam::BlockLduPrecon<Type>::New
(
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict,
    const word fieldName
)
{
    word preconName = getName(dict);

    const dictionary& controls = dict.subOrEmptyDict(preconName+"Coeffs");

    if (matrix.diagonal())
    {
        // No preconditioning for the diagonal matrix
        return autoPtr<BlockLduPrecon<Type>>
        (
            new BlockNoPrecon<Type>
            (
                matrix,
                controls,
                fieldName
            )
        );
    }
    else
    {
        auto ctor =
            ctorTableLookup
            (
                "matrix preconditioner",
                dictionaryConstructorTable_(),
                preconName,
                controls
            );

        return autoPtr<BlockLduPrecon<Type>>
        (
            ctor
            (
                matrix,
                controls,
                fieldName
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::word Foam::BlockLduPrecon<Type>::getName(const dictionary& dict)
{
    return dict.lookup<word>("preconditioner");
}


// ************************************************************************* //
