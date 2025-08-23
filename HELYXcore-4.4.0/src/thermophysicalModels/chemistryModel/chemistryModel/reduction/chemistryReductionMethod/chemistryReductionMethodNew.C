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

#include "chemistryModel/reduction/noChemistryReduction/noChemistryReduction.H"
#include "chemistryModel/reduction/chemistryReductionMethod/chemistryReductionMethod.H"
#include "include/dummyThermo.H"
#include "db/Time/Time.H"
#include "basicThermo/basicThermo.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::autoPtr<Foam::chemistryReductionMethod<ThermoType>>
Foam::chemistryReductionMethod<ThermoType>::New
(
    const IOdictionary& dict,
    chemistryModel<ThermoType>& chemistry
)
{
    if (dict.found("reduction"))
    {
        const dictionary& reductionDict(dict.subDict("reduction"));

        const word methodName(reductionDict.lookup("method"));

        Info<< "Selecting chemistry reduction method " << methodName << endl;

        const word thermoName
        (
            basicThermo::dictName != basicThermo::matDictName
          ? ThermoType::typeName()
          : dummyThermo::typeName()
        );

        const word methodTypeName(methodName + '<' + thermoName + '>');

        // Advanced error message not needed for material properties
        const auto ctor =
            ctorTableLookup
            (
                "chemistryReductionMethodType type",
                dictionaryConstructorTable_(),
                methodTypeName
            );

        autoPtr<chemistryReductionMethod<ThermoType>> crmPtr
        (
            ctor(dict, chemistry)
        );

        chemistry.reduction_ = crmPtr->active();

        return crmPtr;
    }
    else
    {
        return autoPtr<chemistryReductionMethod<ThermoType>>
        (
            new chemistryReductionMethods::none<ThermoType>(dict, chemistry)
        );
    }
}


// ************************************************************************* //
