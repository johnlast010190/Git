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
    (c) 2013-2022 OpenFOAM Foundation
    (c) 2024-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "basicSolidChemistryModel/basicSolidChemistryModel.H"
#include "include/dummyThermo.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicSolidChemistryModel>
Foam::basicSolidChemistryModel::New(const solidMulticomponentThermo& thermo)
{
    IOdictionary chemistryDict
    (
        IOobject
        (
            thermo.phasePropertyName("chemistryProperties"),
            thermo.db().time().constant(),
            thermo.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const dictionary& chemistryTypeDict =
        chemistryDict.subDict("chemistryType");

    Info<< "Selecting chemistry type " << chemistryTypeDict << endl;

    word chemistryTypeName;

    const word solverName
    (
        chemistryTypeDict.found("solver")
      ? chemistryTypeDict.lookup("solver")
      : chemistryTypeDict.found("chemistrySolver")
      ? chemistryTypeDict.lookup("chemistrySolver")
      : chemistryTypeDict.lookup("solver") // error if neither entry is found
    );

    if (basicThermo::dictName != basicThermo::matDictName)
    {
        IOdictionary thermoDict
        (
            IOobject
            (
                basicThermo::dictName,
                thermo.db().time().constant(),
                thermo.db(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        );

        const dictionary& solidThermoTypeDict =
            thermoDict.subDict("thermoType");
        word solidThermoTypeName
        (
            solidThermoTypeDict.lookup<word>("transport") + '<'
          + solidThermoTypeDict.lookup<word>("thermo") + '<'
          + solidThermoTypeDict.lookup<word>("equationOfState") + '<'
          + solidThermoTypeDict.lookup<word>("specie") + ">>,"
          + solidThermoTypeDict.lookup<word>("energy") + ">"
        );

        const dictionary& gasThermoTypeDict =
            thermoDict.subDict("gasThermoType");
        word gasThermoTypeName
        (
            gasThermoTypeDict.lookup<word>("transport") + '<'
          + gasThermoTypeDict.lookup<word>("thermo") + '<'
          + gasThermoTypeDict.lookup<word>("equationOfState") + '<'
          + gasThermoTypeDict.lookup<word>("specie") + ">>,"
          + gasThermoTypeDict.lookup<word>("energy") + ">"
        );

        // Construct the name of the chemistry type from the components
        chemistryTypeName =
            solverName + '<'
          + chemistryTypeDict.lookup<word>("chemistryThermo") + '<'
          + solidThermoTypeName + ',' + gasThermoTypeName + ">>";
    }
    else
    {
        chemistryTypeName =
            solverName + '<'
          + chemistryTypeDict.lookup<word>("chemistryThermo") + '<'
          + dummyThermo::typeName() + ',' + dummyThermo::typeName() + ">>";
    }

    Info<< "chemistryTypeName " << chemistryTypeName << endl;

    thermoConstructorTable::iterator cstrIter =
        thermoConstructorTable_().find(chemistryTypeName);

    if (cstrIter == thermoConstructorTable_().end())
    {
        FatalErrorInFunction
            << "Unknown " << typeName << " type " << nl
            << "chemistryType" << chemistryTypeDict << nl << nl
            << "Valid " << typeName << " types are:"
            << nl << nl;

        // Get the list of all the suitable chemistry packages available
        wordList names;
        for (const auto& p : thermoConstructorTable_())
        {
            names.append(p.first);
        }
        Info<< names << endl;

        // Build a table of the thermo packages constituent parts
        // Note: row-0 contains the names of constituent parts
        List<wordList> validChemistryTypeNameCmpts(names.size() + 1);

        const wordList cmptNames
        ({
            "solver",
            "chemistryThermo",
            "baseChemistry",
            "transport",
            "thermo",
            "equationOfState",
            "specie",
            "energy",
            "transport",
            "thermo",
            "equationOfState",
            "specie",
            "energy"
        });
        const label nCmpts = cmptNames.size();

        validChemistryTypeNameCmpts[0] = cmptNames;

        // Split the thermo package names into their constituent parts
        forAll(names, i)
        {
            validChemistryTypeNameCmpts[i + 1] =
                basicThermo::splitThermoName(names[i], nCmpts);
        }

        // Print the table of available packages
        // in terms of their constituent parts
        printTable(validChemistryTypeNameCmpts, FatalErrorInFunction);

        FatalErrorInFunction << exit(FatalError);
    }

    return autoPtr<basicSolidChemistryModel>(cstrIter->second(thermo));
}


// ************************************************************************* //
