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
    (c) 2015 OpenCFD Ltd.
    (c) 2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "runTimeControl/runTimeCondition/runTimeCondition/runTimeCondition.H"

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::functionObjects::runTimeControls::runTimeCondition>
Foam::functionObjects::runTimeControls::runTimeCondition::New
(
    const word& conditionName,
    const objectRegistry& obr,
    const dictionary& dict,
    stateFunctionObject& state
)
{
    word conditionType(dict.lookup("type"));

    Info<< "Selecting runTimeCondition " << conditionType << endl;

    const auto ctor =
        ctorTableLookup
        (
            "runTimeCondition type",
            dictionaryConstructorTable_(),
            conditionType
        );
    return autoPtr<runTimeCondition>(ctor(conditionName, obr, dict, state));
}


// ************************************************************************* //
