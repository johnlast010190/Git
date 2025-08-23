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
    (c) 2019-2020 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "functionObjects/specieReactionRates/specieReactionRates.H"
#include "fields/volFields/volFields.H"
#include "finiteVolume/fvc/fvcVolumeIntegrate.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(specieReactionRates, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        specieReactionRates,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::specieReactionRates::writeFileHeader
(
    Ostream& os
) const
{
    writeHeader(os, "Specie reaction rates");
    volRegion::writeFileHeader(*this, os);
    writeHeaderValue(os, "nSpecie", chemistryModel_.nSpecie());
    writeHeaderValue(os, "nReaction", chemistryModel_.nReaction());

    writeCommented(os, "Time");
    writeDelimited(os, "Reaction");

    const wordList& speciesNames =
        chemistryModel_.thermo().composition().species();

    forAll(speciesNames, si)
    {
        writeDelimited(os, speciesNames[si] + "_min");
        writeDelimited(os, speciesNames[si] + "_max");
        writeDelimited(os, speciesNames[si] + "_mean");
    }
    writeDelimited(os, "Qdot_min");
    writeDelimited(os, "Qdot_max");
    writeDelimited(os, "Qdot_mean");

    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::specieReactionRates::specieReactionRates
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    volRegion(fvMeshFunctionObject::mesh_, dict),
    writeFile(obr_, name, typeName, dict),
    chemistryModel_
    (
        fvMeshFunctionObject::obr().template lookupObject<basicChemistryModel>
        (
            IOobject::groupName
            (
                "chemistryProperties",
                dict.lookupOrDefault("phaseName", word::null)
            )
        )
    ),
    writeRRFields_
    (
        dict.lookupOrDefault("writeRRFields", List<Tuple2<word, label>>())
    ),
    writeQdotFields_
    (
        dict.lookupOrDefault("writeQdotFields", labelList())
    )
{
    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::specieReactionRates::~specieReactionRates()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::specieReactionRates::read
(
    const dictionary& dict
)
{
    regionFunctionObject::read(dict);

    return true;
}


bool Foam::functionObjects::specieReactionRates::execute()
{
    return true;
}


bool Foam::functionObjects::specieReactionRates::write()
{
    const fvMesh& mesh = this->fvMeshFunctionObject::mesh_;
    const label nSpecie = chemistryModel_.nSpecie();
    const label nReaction = chemistryModel_.nReaction();

    // Region volume
    const scalar V = this->V();

    for (label reactioni=0; reactioni<nReaction; reactioni++)
    {
        writeTime(file());
        file() << token::TAB << reactioni;

        // Sum heat release for this reaction
        volScalarField::Internal Qdot
        (
            IOobject
            (
                "Qdot",
                mesh.time().timeName(),
                obr()
            ),
            mesh,
            dimensionedScalar(dimPower/dimVolume, 0)
        );

        const PtrList<volScalarField::Internal> RR
        (
            chemistryModel_.reactionRR(reactioni)
        );

        for (label speciei=0; speciei<nSpecie; speciei++)
        {
            word speciesName =
                chemistryModel_.thermo().composition().species()[speciei];

            scalar sumVRRi = 0;
            scalar minVRRi = GREAT;
            scalar maxVRRi = -GREAT;

            if (isNull(cellIDs()))
            {
                sumVRRi = fvc::domainIntegrate(RR[speciei]).value();
                minVRRi = gMin(RR[speciei]);
                maxVRRi = gMax(RR[speciei]);
            }
            else
            {
                sumVRRi = gSum(scalarField(mesh.V()*RR[speciei], cellIDs()));
                minVRRi = gMin(scalarField(RR[speciei], cellIDs()));
                maxVRRi = gMax(scalarField(RR[speciei], cellIDs()));
            }

            file() << token::TAB << minVRRi;
            file() << token::TAB << maxVRRi;
            file() << token::TAB << sumVRRi/V;

            forAll(writeRRFields_, i)
            {
                if
                (
                    writeRRFields_[i]
                 == Tuple2<word, label>(speciesName, reactioni)
                )
                {
                    OStringStream name;
                    name << "RR_" << speciesName << "_" << reactioni;
                    if (!obr().template foundObject<volScalarField>(name.str()))
                    {
                        regIOobject::store
                        (
                            new volScalarField
                            (
                                IOobject
                                (
                                    name.str(),
                                    mesh.time().timeName(),
                                    obr(),
                                    IOobject::NO_READ,
                                    IOobject::AUTO_WRITE
                                ),
                                mesh,
                                RR[speciei].dimensions(),
                                zeroGradientFvPatchField<scalar>::typeName
                            )
                        );
                    }
                    volScalarField& RRField =
                        obr().template lookupObjectRef<volScalarField>
                        (
                            name.str()
                        );
                    RRField.ref() = RR[speciei];
                    RRField.correctBoundaryConditions();
                }
            }
        }

        // Write heat release
        scalar sumVQdoti = 0;
        scalar minVQdoti = GREAT;
        scalar maxVQdoti = -GREAT;

        if (isNull(cellIDs()))
        {
            minVQdoti = gMin(Qdot);
            maxVQdoti = gMax(Qdot);
            sumVQdoti = fvc::domainIntegrate(Qdot).value();
        }
        else
        {
            sumVQdoti = gSum
            (
                scalarField(fvMeshFunctionObject::mesh_.V()*Qdot, cellIDs())
            );
            minVQdoti = gMin(scalarField(Qdot, cellIDs()));
            maxVQdoti = gMax(scalarField(Qdot, cellIDs()));
        }

        file() << token::TAB << minVQdoti;
        file() << token::TAB << maxVQdoti;
        file() << token::TAB << sumVQdoti/V;

        forAll(writeQdotFields_, i)
        {
            if (writeQdotFields_[i] == reactioni)
            {
                OStringStream name;
                name << "Qdot_" << reactioni;
                if (!obr().template foundObject<volScalarField>(name.str()))
                {
                    regIOobject::store
                    (
                        new volScalarField
                        (
                            IOobject
                            (
                                name.str(),
                                mesh.time().timeName(),
                                obr(),
                                IOobject::NO_READ,
                                IOobject::AUTO_WRITE
                            ),
                            mesh,
                            Qdot.dimensions(),
                            zeroGradientFvPatchField<scalar>::typeName
                        )
                    );
                }
                volScalarField& QdotField =
                    obr().template lookupObjectRef<volScalarField>(name.str());
                QdotField.ref() = Qdot;
                QdotField.correctBoundaryConditions();
            }
        }

        file() << nl;
    }

    file() << endl;

    return true;
}


// ************************************************************************* //
