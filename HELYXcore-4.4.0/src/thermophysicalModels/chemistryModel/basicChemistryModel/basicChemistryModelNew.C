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
    (c) 2012-2022 OpenFOAM Foundation
    (c) 2021-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "basicChemistryModel/basicChemistryModel.H"
#include "include/dummyThermo.H"
#include "basicThermo/basicThermo.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicChemistryModel> Foam::basicChemistryModel::New
(
    const fluidMulticomponentThermo& thermo
)
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

    const word solverName
    (
        chemistryTypeDict.found("solver")
      ? chemistryTypeDict.lookup("solver")
      : chemistryTypeDict.found("chemistrySolver")
      ? chemistryTypeDict.lookup("chemistrySolver")
      : chemistryTypeDict.lookup("solver") // error if neither entry is found
    );

    const word methodName
    (
        chemistryTypeDict.lookupOrDefault<word>("method", "chemistryModel")
    );

    dictionary chemistryTypeDictNew;
    chemistryTypeDictNew.add("solver", solverName);
    chemistryTypeDictNew.add("method", methodName);

    Info<< "Selecting chemistry solver " << chemistryTypeDictNew << endl;

    const word thermoName
    (
        basicThermo::dictName != basicThermo::matDictName
      ? thermo.thermoName()
      : dummyThermo::typeName()
    );

    const word chemSolverNameName
    (
        solverName + '<' + methodName + '<' + thermoName + ">>"
    );

    typename thermoConstructorTable::iterator cstrIter =
        thermoConstructorTable_().find(chemSolverNameName);

    if (cstrIter == thermoConstructorTable_().end())
    {
        FatalErrorInFunction
            << "Unknown " << typeName_() << " type " << solverName << endl
            << endl;

        wordList names;
        for (const auto& p : thermoConstructorTable_())
        {
            names.append(p.first);
        }

        wordList thisCmpts({word::null, word::null});
        thisCmpts.append(basicThermo::splitThermoName(thermoName, 5));

        List<wordList> validNames;
        validNames.append(wordList({"solver", "method"}));
        forAll(names, i)
        {
            const wordList cmpts(basicThermo::splitThermoName(names[i], 7));

            if (SubList<word>(cmpts, 5, 2) == SubList<word>(thisCmpts, 5, 2))
            {
                validNames.append(SubList<word>(cmpts, 2));
            }
        }

        FatalErrorInFunction
            << "Valid " << validNames[0][0] << '/' << validNames[0][1]
            << "combinations for this thermodynamic model are:"
            << endl << endl;
        printTable(validNames, FatalErrorInFunction);

        FatalErrorInFunction << endl;

        List<wordList> validCmpts;
        validCmpts.append
        (
            wordList
            ({
                "solver",
                "method",
                "transport",
                "thermo",
                "equationOfState",
                "specie",
                "energy"
            })
        );
        forAll(names, i)
        {
            validCmpts.append(basicThermo::splitThermoName(names[i], 7));
        }

        FatalErrorInFunction
            << "All " << validCmpts[0][0] << '/' << validCmpts[0][1]
            << "/thermodynamics combinations are:"
            << endl << endl;
        printTable(validCmpts, FatalErrorInFunction);

        FatalErrorInFunction << exit(FatalError);
    }

    return cstrIter->second(thermo);
}

// ************************************************************************* //
