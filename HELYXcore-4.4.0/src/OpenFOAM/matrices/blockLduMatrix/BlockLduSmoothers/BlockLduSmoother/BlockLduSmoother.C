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

Class
    BlockLduSmoother

Description
    Block LDU matrix smoother virtual base class

\*---------------------------------------------------------------------------*/

#include "matrices/blockLduMatrix/BlockLduSmoothers/BlockLduSmoother/BlockLduSmoother.H"

template<class Type>
class BlockNoSmoother;

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::BlockLduSmoother<Type>> Foam::BlockLduSmoother<Type>::New
(
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict
)
{
    word smootherName;

    // Handle primitive or dictionary entry
    const entry& e = dict.lookupEntry("smoother", false, false);
    if (e.isDict())
    {
        smootherName = e.dict().lookup<word>("smoother");
    }
    else
    {
        e.stream() >> smootherName;
    }

    // Not (yet?) needed:
    // const dictionary& controls = e.isDict() ? e.dict() : dictionary::null;

    auto ctor =
        ctorTableLookup
        (
            "matrix smoother",
            dictionaryConstructorTable_(),
            smootherName,
            dict
        );

    return autoPtr<BlockLduSmoother<Type>>
    (
        ctor
        (
            matrix,
            dict
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::word Foam::BlockLduSmoother<Type>::getName(const dictionary& dict)
{
    word name;

    // handle primitive or dictionary entry
    const entry& e = dict.lookupEntry("preconditioner", false, false);
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


// ************************************************************************* //
