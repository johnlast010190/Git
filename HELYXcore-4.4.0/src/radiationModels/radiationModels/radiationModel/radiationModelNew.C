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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "radiationModels/radiationModel/radiationModel.H"
#include "fields/volFields/volFields.H"
#include "db/IOobjectList/IOobjectList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::radiationModel>
Foam::radiationModel::New
(
    const volScalarField& T, const dimensionedScalar& TRef
)
{
    // check if radiation model already exists
    if (T.mesh().foundObject<IOdictionary>(radiationModel::typeName))
    {
        FatalErrorInFunction
            << "Radiation model already exists in registry, please check your setup!"
            << nl << "You have probably added a radiation model via fvOption"
            << " to a solver that already uses a radiation model."
            << exit(FatalError);
    }

    IOobject radIO
    (
        dictName,
        T.db().time().constant(),
        T.db(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        false
    );

    word modelType("none");
    if (radIO.typeHeaderOk<IOdictionary>(true))
    {
        modelType = IOdictionary(radIO).lookup<word>("radiationModel");
    }
    else
    {
        Info<< "Radiation model not active: radiationProperties not found"
            << endl;
    }

    Info<< "Selecting radiationModel " << modelType << endl;

    const auto ctor =
        ctorTableLookup
        (
            "radiationModel type",
            TConstructorTable_(),
            modelType
        );
    return autoPtr<radiationModel>(ctor(T, TRef));
}


Foam::autoPtr<Foam::radiationModel>
Foam::radiationModel::New
(
    const dictionary& dict,
    const volScalarField& T,
    const dimensionedScalar& TRef
)
{
    const word modelType(dict.lookup("radiationModel"));

    Info<< "Selecting radiationModel " << modelType << endl;

    const auto ctor =
        ctorTableLookup
        (
            "radiation model",
            dictionaryConstructorTable_(),
            modelType
        );
    return autoPtr<radiationModel>(ctor(dict, T, TRef));
}


Foam::radiationModel& Foam::radiationModel::lookupOrCreate
(
    const volScalarField& T,
    const dimensionedScalar& TRef
)
{
    if (!T.db().foundObject<radiationModel>(dictName))
    {
        if (radiationModel::debug)
        {
            InfoInFunction
                << "constructing radiation model " << dictName
                << " for region " << T.db().name() << endl;
        }
        regIOobject::store(New(T, TRef).ptr());
    }

    return T.db().lookupObjectRef<radiationModel>(dictName);
}


// ************************************************************************* //
