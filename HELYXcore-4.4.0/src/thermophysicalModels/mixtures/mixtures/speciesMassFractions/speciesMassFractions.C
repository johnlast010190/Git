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
    (c) 2011-2020 OpenFOAM Foundation
    (c) 2021-2025 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "speciesMassFractions.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(speciesMassFractions, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void Foam::speciesMassFractions::createFields
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const objectRegistry& obr,
    const word& phaseName
)
{
    word YdefaultName(IOobject::groupName("Ydefault", phaseName));
    tmp<volScalarField> tYdefault;
    const fvMesh& mesh = fvSolutionRegistry::getMesh(obr);
    forAll(species_, i)
    {
        IOobject header
        (
            IOobject::groupName(species_[i], phaseName),
            obr.time().timeName(),
            obr,
            IOobject::NO_READ
        );

        // Check if field exists and can be looked up or read
        word YName(IOobject::groupName(species_[i], phaseName));
        if
        (
            obr.foundObject<volScalarField>(YName)
         && obr.lookupObject<volScalarField>(YName).ownedByRegistry()
        )
        {
            // Take ownership
            volScalarField& Y = obr.lookupObjectRef<volScalarField>(YName);
            Y.release();
            Y_.set(i, &Y);
        }
        else if (header.typeHeaderOk<volScalarField>(true))
        {
            DebugInformation
                << "speciesMassFractions: reading " << species_[i]
                << endl;

            Y_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        YName,
                        obr.time().timeName(),
                        obr,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        }
        else
        {
            // Read Ydefault if not already read
            if (!tYdefault.valid())
            {
                if (!obr.foundObject<volScalarField>(YdefaultName))
                {
                    IOobject timeIO
                    (
                        YdefaultName,
                        obr.time().timeName(),
                        obr,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    );

                    IOobject constantIO
                    (
                        YdefaultName,
                        obr.time().constant(),
                        obr,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    );

                    IOobject time0IO
                    (
                        YdefaultName,
                        Time::timeName(0),
                        obr,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    );

                    if (timeIO.typeHeaderOk<volScalarField>(true))
                    {
                        tYdefault = new volScalarField(timeIO, mesh);
                    }
                    else if (constantIO.typeHeaderOk<volScalarField>(true))
                    {
                        tYdefault = new volScalarField(constantIO, mesh);
                    }
                    else
                    {
                        tYdefault = new volScalarField(time0IO, mesh);
                    }
                }
            }

            Y_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName(species_[i], phaseName),
                        obr.time().timeName(),
                        obr,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    obr.lookupObject<volScalarField>(YdefaultName)
                )
            );
        }
    }

    correctMassFractions();
}


void Foam::speciesMassFractions::setDefaultSpecie(const dictionary& dict)
{
    if (defaultSpecie_ == "none")
    {
        defaultSpecie_ = word::null;
    }
    else if (defaultSpecie_ == word::null)
    {
        // Since the estimate of the inertSpecie doesn't have to be
        // exact the simple sum over field should be enough to compare
        // the fields.
        scalarField volSumY(Y_.size());
        forAll(Y_, i)
        {
            volSumY[i] =
                gSum(Y_[i].primitiveField()) + gSum(Y_[i].boundaryField());
        }
        defaultSpecie_ = Y_[findMax(volSumY)].member();
    }
    if (defaultSpecie_ != word::null)
    {
        if (!species().found(defaultSpecie_))
        {
            FatalIOErrorInFunction(dict)
                << "Inert/default species " << defaultSpecie_
                << " not found in available species " << species()
                << exit(FatalIOError);
        }
        else
        {
            defaultSpecieIndex_ = species()[defaultSpecie_];
        }
    }
}


Foam::speciesMassFractions::speciesMassFractions
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const objectRegistry& obr,
    const word& phaseName
)
:
    volumeMassFractions(obr, phaseName),
    species_(specieNames),
    defaultSpecie_
    (
        thermoDict.optionalSubDict(phaseName).lookupOrDefault<word>
        (
            wordList({"inertSpecie", "defaultSpecie"}),
            word::null
        )
    ),
    defaultSpecieIndex_(-1),
    active_(species_.size(), true),
    Y_(species_.size())
{
    if (species_.size())
    {
        createFields(thermoDict, specieNames, obr, phaseName);
        setDefaultSpecie(thermoDict);
    }
}


Foam::speciesMassFractions::speciesMassFractions
(
    const dictionary& thermoDict,
    const objectRegistry& obr,
    const word& phaseName
)
:
    volumeMassFractions(obr, phaseName),
    species_
    (
        (phaseName != word::null)
      ? thermoDict.optionalSubDict
        (
            phaseName
        ).lookupOrDefault<wordList>("species", {})
      : thermoDict.lookupOrDefault<wordList>("species", {})
    ),
    defaultSpecie_
    (
        thermoDict.optionalSubDict(phaseName).lookupOrDefault<word>
        (
            wordList({"inertSpecie", "defaultSpecie"}),
            word::null
        )
    ),
    defaultSpecieIndex_(-1),
    active_(species_.size(), true),
    Y_(species_.size())
{
    if (species_.size())
    {
        createFields(thermoDict, species_, obr, phaseName);
        setDefaultSpecie(thermoDict);
    }
}


void Foam::speciesMassFractions::correctMassFractions()
{
    if (species_.size())
    {
        // Multiplication by 1.0 changes Yt patches to "calculated"
        volScalarField Yt("Yt", 1.0*Y_[0]);

        for (label i = 1; i < Y_.size(); i++)
        {
            Yt += Y_[i];
        }

        if (min(Yt).value() < ROOTVSMALL)
        {
            FatalErrorInFunction
                << "Sum of species mass fractions is zero or less";
            if (min(Yt.primitiveField()) < ROOTVSMALL)
            {
                FatalErrorInFunction
                    << " in internal field";
            }
            else
            {
                forAll(Yt.boundaryField(), patchi)
                {
                    if (min(Yt.boundaryField()[patchi]) < ROOTVSMALL)
                    {
                        FatalErrorInFunction
                            << " in patch field "
                            << Yt.boundaryField()[patchi].patch().name();
                    }
                }
            }
            FatalErrorInFunction
                << " (fields:";
            forAll(Y_, Yi)
            {
                FatalErrorInFunction << " " << Y_[Yi].name();
            }
            FatalErrorInFunction << ")." << nl << exit(FatalError);
        }

        forAll(Y_, n)
        {
            Y_[n] /= Yt;
        }
    }
}


void Foam::speciesMassFractions::normalise()
{
    if (species_.size())
    {
        volScalarField Yt(0.0*Y_[0]);
        forAll(Y_, i)
        {
            if (solve(i))
            {
                Y_[i].max(0.0);
                Yt += Y_[i];
            }
        }

        Y_[defaultSpecieIndex_] = scalar(1) - Yt;
        Y_[defaultSpecieIndex_].max(0.0);
    }
}

// ************************************************************************* //
