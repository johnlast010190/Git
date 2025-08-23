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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "submodels/Kinematic/InjectionModel/InjectionModel/InjectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::autoPtr<Foam::InjectionModel<CloudType>>
Foam::InjectionModel<CloudType>::New
(
    const dictionary& dict,
    CloudType& owner
)
{
    const word modelType(dict.lookup("injectionModel"));

    Info<< "Selecting injection model " << modelType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTable_().find(modelType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown injection model type "
            << modelType << nl << nl
            << "Valid injection model types are:" << nl
            << dictionaryConstructorTable_().sortedToc() << exit(FatalError);
    }

    return autoPtr<InjectionModel<CloudType>>(cstrIter->second(dict, owner));
}


template<class CloudType>
Foam::autoPtr<Foam::InjectionModel<CloudType>>
Foam::InjectionModel<CloudType>::New
(
    const dictionary& dict,
    const word& modelName,
    const word& modelType,
    CloudType& owner
)
{
    Info<< "Selecting injection model " << modelType << endl;

    const auto ctor =
        ctorTableLookup
        (
            "injection model type",
            dictionaryConstructorTable_(),
            modelType
        );
    return
        autoPtr<InjectionModel<CloudType>>
        (
            ctor
            (
                dict,
                owner,
                modelName
            )
        );
}


// ************************************************************************* //
