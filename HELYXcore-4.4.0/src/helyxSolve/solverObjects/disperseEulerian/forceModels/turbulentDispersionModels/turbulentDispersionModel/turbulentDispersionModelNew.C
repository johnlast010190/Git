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

#include "turbulentDispersionModel.H"
#include "solverObjects/disperseEulerian/phase/phase.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::decoupledEulerian::turbulentDispersionModel>
Foam::decoupledEulerian::turbulentDispersionModel::New
(
    const dictionary& dict,
    const phase& phase,
    const bool registerObject
)
{
    word turbulentDispersionModelType(dict.lookup("turbulentDispersionModel"));

    Info<< "Selecting turbulentDispersionModel for "
        << phase.name() << ": " << turbulentDispersionModelType << endl;

    const auto ctor =
        ctorTableLookup
        (
            "turbulentDispersionModelType type",
            dictionaryConstructorTable_(),
            turbulentDispersionModelType
        );
    return ctor(dict, phase, registerObject);
}


// ************************************************************************* //
