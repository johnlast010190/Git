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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "matrices/lduMatrix/solvers/GAMG/interfaces/GAMGInterface/GAMGInterface.H"
#include "matrices/lduMatrix/solvers/GAMG/GAMGAgglomerations/GAMGAgglomeration/GAMGAgglomeration.H"
#include "matrices/lduMatrix/lduMatrix/lduMatrix.H"


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::GAMGInterface> Foam::GAMGInterface::New
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    const lduInterface& fineInterface,
    const labelField& localRestrictAddressing,
    const labelField& neighbourRestrictAddressing,
    const label fineLevelIndex,
    const label coarseComm
)
{
    const word coupleType(fineInterface.type());

    const auto ctor =
        ctorTableLookup
        (
            "GAMGInterface type",
            lduInterfaceConstructorTable_(),
            coupleType
        );
    return autoPtr<GAMGInterface>
    (
        ctor
        (
            index,
            coarseInterfaces,
            fineInterface,
            localRestrictAddressing,
            neighbourRestrictAddressing,
            fineLevelIndex,
            coarseComm
        )
    );
}


Foam::autoPtr<Foam::GAMGInterface> Foam::GAMGInterface::New
(
    const word& coupleType,
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    Istream& is
)
{
    const auto ctor =
        ctorTableLookup
        (
            "GAMGInterface type",
            IstreamConstructorTable_(),
            coupleType
        );
    return autoPtr<GAMGInterface>(ctor(index, coarseInterfaces, is));
}


// ************************************************************************* //
