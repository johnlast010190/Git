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
    (c) 2016-2022 OpenFOAM Foundation
    (c) 2021-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "chemistryModel/tabulation/chemistryTabulationMethod/chemistryTabulationMethod.H"
#include "include/dummyThermo.H"
#include "db/Time/Time.H"
#include "basicThermo/basicThermo.H"
#include "chemistryModel/tabulation/noChemistryTabulation/noChemistryTabulation.H"


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::chemistryTabulationMethod>
Foam::chemistryTabulationMethod::New
(
    const IOdictionary& dict,
    const odeChemistryModel& chemistry
)
{
    if (dict.found("tabulation"))
    {
        const dictionary& tabulationDict = dict.subDict("tabulation");

        const word methodName(tabulationDict.lookup("method"));

        Info<< "Selecting chemistry tabulation method " << methodName << endl;

        const auto ctor =
            ctorTableLookup
            (
                "chemistryTabulationMethodType type",
                dictionaryConstructorTable_(),
                methodName
            );
        return autoPtr<chemistryTabulationMethod>(ctor(dict, chemistry));
    }
    else
    {
        return autoPtr<chemistryTabulationMethod>
        (
            new chemistryTabulationMethods::none(dict, chemistry)
        );
    }
}


// ************************************************************************* //
