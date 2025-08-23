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
    (c) 2012-2019 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "interpolations/interpolationWeights/interpolationWeights/interpolationWeights.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(interpolationWeights, 0);
defineRunTimeSelectionTable(interpolationWeights, word);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

interpolationWeights::interpolationWeights
(
    const scalarField& samples
)
:
    samples_(samples)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<interpolationWeights> interpolationWeights::New
(
    const word& type,
    const scalarField& samples
)
{
    if (debug)
    {
        InfoInFunction
            << "Selecting interpolationWeights "
            << type << endl;
    }

    const auto ctor =
        ctorTableLookup
        (
            "interpolationWeights type",
            wordConstructorTable_(),
            type
        );
    return autoPtr<interpolationWeights>(ctor(samples));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

interpolationWeights::~interpolationWeights()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
