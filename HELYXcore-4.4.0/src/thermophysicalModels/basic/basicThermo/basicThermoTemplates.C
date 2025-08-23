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
    (c) 2012-2017 OpenFOAM Foundation
    (c) 2020-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "basicThermo/basicThermo.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Thermo, class Table>
typename Table::iterator Foam::basicThermo::lookupThermo
(
    const dictionary& thermoTypeDict,
    Table& tablePtr,
    const int nCmpt,
    const char* cmptNames[],
    const word& thermoTypeName,
    const word& phaseName
)
{
    // Lookup the thermo package
    typename Table::iterator cstrIter = tablePtr.find(thermoTypeName);

    // Print error message if package not found in the table
    if (cstrIter == tablePtr.end())
    {
        FatalErrorInFunction
            << "Unknown " << Thermo::typeName << " type " << nl
            << "thermoType" << thermoTypeDict << nl << nl
            << "Valid " << Thermo::typeName << " types are:"
            << nl << nl;

        // Get the list of all the suitable thermo packages available
        wordList validThermoTypeNames;
        for (const auto& p : tablePtr)
        {
            validThermoTypeNames.append(p.first);
        }

        // Build a table of the thermo packages constituent parts
        // Note: row-0 contains the names of constituent parts
        DynamicList<wordList> validThermoTypeNameCmpts;;

        // Set row zero to the column headers
        validThermoTypeNameCmpts.append(wordList(nCmpt));
        forAll(validThermoTypeNameCmpts[0], i)
        {
            validThermoTypeNameCmpts[0][i] = cmptNames[i];
        }

        // Split the thermo package names into their constituent parts and add
        // them to the table, removing any incompatible entries from the list
        forAll(validThermoTypeNames, i)
        {
            const wordList names
            (
                Thermo::splitThermoName(validThermoTypeNames[i], nCmpt)
            );

            if (names.size())
            {
                validThermoTypeNameCmpts.append(names);
            }
        }

        // Print the table of available packages
        printTable(validThermoTypeNameCmpts, FatalError);

        FatalError<< exit(FatalError);
    }

    return cstrIter;
}


template<class Thermo, class Table>
typename Table::iterator Foam::basicThermo::lookupThermo
(
    const dictionary& thermoDict,
    Table& tablePtr,
    const word& phaseName
)
{
    if (thermoDict.isDict("thermoType"))
    {
        const dictionary& thermoTypeDict = thermoDict.subDict("thermoType");

        Info<< "Selecting thermodynamics package " << thermoTypeDict << endl;

        const int nCmpt = 7;
        const char* cmptNames[nCmpt] =
        {
            "type",
            "mixture",
            "transport",
            "thermo",
            "equationOfState",
            "specie",
            "energy"
        };

        // Construct the name of the thermo package from the components
        const word thermoTypeName
        (
            thermoTypeDict.lookup<word>("type") + '<'
          + thermoTypeDict.lookup<word>("mixture") + '<'
          + thermoTypeDict.lookup<word>("transport") + '<'
          + thermoTypeDict.lookup<word>("thermo") + '<'
          + thermoTypeDict.lookup<word>("equationOfState") + '<'
          + thermoTypeDict.lookup<word>("specie") + ">>,"
          + thermoTypeDict.lookup<word>("energy") + ">>>"
        );

        return lookupThermo<Thermo, Table>
        (
            thermoTypeDict,
            tablePtr,
            nCmpt,
            cmptNames,
            thermoTypeName
        );
    }
    else
    {
        if (thermoDict.found("thermoType"))
        {
            FatalErrorInFunction
                << "thermoType specification isn't supported option." << nl
                << exit(FatalError);
        }

        word materialType_
        (
            thermoDict.optionalSubDict(phaseName).lookup("materialType")
        );
        word extraPhaseMessage;
        word phaseDashes;
        if (phaseName != word::null)
        {
            extraPhaseMessage = " \"" + phaseName + "\"";
            phaseDashes = std::string(phaseName.size() + 3, '-');
        }
        Info<< nl
            << "Material settings" << extraPhaseMessage << nl
            << "-----------------" << phaseDashes << nl << nl
            << "Selecting materials package \"" << materialType_
            << "\"." << endl;

        //TODO: In future even without multiple species this could be reacting
        // mixture. That way no "hack" is required to load up mat. props
        if
        (
            materialType_ == "fluid"
         && thermoDict.optionalSubDict(phaseName).found("species")
        )
        {
            materialType_ = "reactingFluid";
        }
        else if
        (
            materialType_ == "solid"
         && thermoDict.optionalSubDict(phaseName).found("species")
        )
        {
            materialType_ = "reactingSolid";
        }

        typename Table::iterator cstrIter = tablePtr.find(materialType_);

        if (cstrIter == tablePtr.end())
        {
            FatalErrorInFunction
                << "Unknown " << Thermo::typeName << " materialType "
                << materialType_ << nl << nl
                << exit(FatalError);
        }

        return cstrIter;
    }
}


template<class Thermo>
Foam::autoPtr<Thermo> Foam::basicThermo::New
(
    const objectRegistry& obr,
    const word& phaseName
)
{
    IOdictionary thermoDict(loadDictionary(obr, phaseName, false));

    typename Thermo::objectRegistryConstructorTable::iterator cstrIter =
        lookupThermo<Thermo, typename Thermo::objectRegistryConstructorTable>
        (
            thermoDict,
            Thermo::objectRegistryConstructorTable_(),
            phaseName
        );

    return autoPtr<Thermo>(cstrIter->second(obr, phaseName));
}


template<class Thermo>
Foam::autoPtr<Thermo> Foam::basicThermo::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName
)
{
    typename Thermo::dictionaryConstructorTable::iterator cstrIter =
        lookupThermo<Thermo, typename Thermo::dictionaryConstructorTable>
        (
            dict,
            Thermo::dictionaryConstructorTable_(),
            phaseName
        );

    return autoPtr<Thermo>(cstrIter->second(obr, dict, phaseName));
}


// ************************************************************************* //
