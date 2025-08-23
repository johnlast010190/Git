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
    (c) 2016-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "multiphaseMulticomponentState.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "phasePairKey/phasePairKey.H"
#include "multiphaseThermo/multiphaseThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace stateFunctions
{

defineTypeNameAndDebug(multiphaseMulticomponentState, 0);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void multiphaseMulticomponentState::initialiseMaterials()
{
    dictionary materialProperties;

    dictionary inDict = input().subOrEmptyDict("materialProperties");

    // To avoid confusion, complain if user tried to directly specify a
    // system/materialProperties file
    if (input().subOrEmptyDict("system").isDict("materialProperties"))
    {
        IOWarningInFunction
        (
            input().subOrEmptyDict("system").subDict("materialProperties")
        )
            << "'materialProperties' dictionary was specified inside 'system'. "
            << "It is advisable to use the 'materialProperties' section "
            << "(one level up) instead."
            << nl << endl;
    }

    // First merge any user-specified entries
    materialProperties.merge(inDict, true);

    // Merge and overwrite from defaults, then re-merge from input
    materialProperties.merge
    (
        defaults().subDict("materialProperties").subDict("referenceValues")
    );
    if (inDict.found("referenceFields"))
    {
        materialProperties.subDict("referenceFields").merge
        (
            inDict.subOrEmptyDict("referenceFields")
        );
    }

    if (defaults().subDict("materialProperties").found("buoyancyModel"))
    {
        materialProperties.add
        (
            "buoyancyModel",
            defaults().subDict("materialProperties").lookup("buoyancyModel")
        );
        materialProperties.add("buoyancyModelCoeffs", dictionary(), false);
        materialProperties.subDict("buoyancyModelCoeffs").merge
        (
            defaults().subDict("materialProperties").subOrEmptyDict
            (
                "buoyancyModelCoeffs"
            )
        );
        materialProperties.subDict("buoyancyModelCoeffs").merge
        (
            inDict.subOrEmptyDict
            (
                "buoyancyModelCoeffs"
            )
        );
    }

    wordList phases(phaseNames());
    if (phases.size() > 1)
    {
        // Multiphase case

        if
        (
            inDict.lookupOrDefault("materialType", word("multiphase"))
         != "multiphase"
        )
        {
            FatalIOErrorInFunction(inDict)
                << "'multiphase' was expected as material type"
                << nl << exit(FatalError);
        }

        // Add required entries if not already present
        materialProperties.add("phases", phases);
        materialProperties.add("materialType", "multiphase");
        materialProperties.add("mixture", "volumeMixture");
        materialProperties.add("energy", "sensibleInternalEnergy");

        forAll(phases, phasei)
        {
            materialProperties.add
            (
                phases[phasei], dictionary(), false
            );

            initialisePhase
            (
                phases[phasei],
                materialProperties.subDict(phases[phasei]),
                inDict.subOrEmptyDict(phases[phasei])
            );
        }
    }
    else
    {
        // Single phase case

        if
        (
            inDict.lookupOrDefault("materialType", word("fluid"))
         == "multiphase"
        )
        {
            FatalIOErrorInFunction(inDict)
                << "'multiphase' was not expected as material type"
                << nl << exit(FatalError);
        }

        // Add required entries if not already present
        materialProperties.add("energy", "sensibleEnthalpy");

        initialisePhase(word::null, materialProperties, inDict);
    }

    system().add("materialProperties", dictionary(), false);
    system().subDict("materialProperties").merge(materialProperties);

}


void multiphaseMulticomponentState::initialisePhase
(
    const word& phaseName,
    dictionary& materialProperties,
    const dictionary& inDict
)
{
    // Called from initialiseMaterials() for each phase or for only phase

    wordList componentNames =
        phaseName == word::null
      ? materialNames()
      : inDict.lookupOrDefault("species", wordList(1, phaseName));

    if (componentNames.size() > 1)
    {
        // Multi-component case

        materialProperties.add("species", componentNames, true);

        materialProperties.add("mixture", "standardMixture");

        forAll(componentNames, specI)
        {
            materialProperties.add
            (
                componentNames[specI], dictionary(), false
            );

            initialiseComponent
            (
                componentNames[specI],
                materialProperties.subDict(componentNames[specI]),
                inDict.subOrEmptyDict(componentNames[specI]),
                inDict
            );
        }

        if (inDict.found("inertSpecie"))
        {
            materialProperties.add
            (
                "inertSpecie",
                inDict.lookup<word>("inertSpecie")
            );
        }
    }
    else
    {
        if (phaseName == word::null)
        {
            // Single-phase/single component case
            // Remove existing subdict that may have been merged in from user
            // input, and re-merge into the parent dict
            materialProperties.remove(componentNames[0]);
            materialProperties.merge(inDict.subOrEmptyDict(componentNames[0]));

            initialiseComponent
            (
                componentNames[0],
                materialProperties,
                inDict.subOrEmptyDict(componentNames[0]),
                inDict
            );
        }
        else
        {
            // Multi-phase/single component case
            initialiseComponent
            (
                componentNames[0],
                materialProperties,
                inDict,
                dictionary()
            );
        }
    }
}


void multiphaseMulticomponentState::initialiseComponent
(
    const word& componentName,
    dictionary& materialProperties,
    const dictionary& inDict,
    const dictionary& inDictParent
)
{
    // Called from initialisePhase for each component or for only component
    materialProperties.merge(inDict);
    // Defaults may contain noise, so don't merge directly
    dictionary allInput;
    allInput.merge
    (
        defaults().subDict("materialProperties").subOrEmptyDict
        (
            componentName + "Material"
        )
    );
    allInput.merge(inDict);
    word matType;
    if (!allInput.found("materialType") && inDictParent.found("materialType"))
    {
        matType = inDictParent.lookup<word>("materialType");
    }
    else
    {
        matType = allInput.lookup<word>("materialType");
        materialProperties.add("materialType", matType);
    }
    materialProperties.merge(speciesProperties(matType, allInput));
}


dictionary multiphaseMulticomponentState::speciesProperties
(
    const word& matType,
    const dictionary& matSpecies
)
{
    dictionary species;
    wordList properties
    ({
        "equationOfState",
        "thermodynamics",
        "kappaModel"
    });
    wordList dictsToAdd
    ({
        "equationOfStateCoeffs",
        "thermodynamicsCoeffs",
        "kappaModelCoeffs"
    });

    if (matType != "solid")
    {
        properties.append("muModel");
        dictsToAdd.append("muModelCoeffs");
    }

    if (matSpecies.found("elements"))
    {
        species.add("elements", matSpecies.subDict("elements"));
    }
    forAll(properties, i)
    {
        species.add(properties[i], matSpecies.lookup<word>(properties[i]));
    }
    if (matSpecies.found("buoyancyModel"))
    {
        species.add("buoyancyModel", matSpecies.lookup<word>("buoyancyModel"));
        if (matSpecies.found("buoyancyModelCoeffs"))
        {
            species.add
            (
                "buoyancyModelCoeffs",
                matSpecies.subDict("buoyancyModelCoeffs")
            );
        }
    }
    species.add("molWeight", matSpecies.lookup<scalar>("molWeight"));

    // Allow to redefine default R for a species
    if (matSpecies.found("R"))
    {
        species.add("R", matSpecies.lookup<scalar>("R"));
    }

    forAll(dictsToAdd, i)
    {
        species.add(dictsToAdd[i], matSpecies.subOrEmptyDict(dictsToAdd[i]));
    }

    return species;
}


void multiphaseMulticomponentState::addAlphaTypeGrad(const word& phaseName)
{
    if (input().found("system"))
    {
        if (input().subDict("system").found("fvSchemes"))
        {
            dictionary& fvSchemes
            (
                input().subDict("system").subDict("fvSchemes")
            );
            if (fvSchemes.found("gradSchemes"))
            {
                if (fvSchemes.subDict("gradSchemes").found("grad(alpha)"))
                {
                    system().add(word("fvSchemes"), dictionary(), false);
                    system().subDict("fvSchemes")
                        .add(word("gradSchemes"), dictionary(), false);

                    system().subDict("fvSchemes").subDict("gradSchemes").add
                    (
                        word("grad(alpha."+phaseName +")"),
                        fvSchemes.subDict("gradSchemes").lookup("grad(alpha)"),
                        false
                    );
                }
            }
        }
    }
}


void multiphaseMulticomponentState::mergePhaseModels(const word& modelName)
{
    dictionary& matDict = system().subDict("materialProperties");
    wordList phases(phaseNames());

    forAll(phases, pi)
    {
        dictionary& dict = matDict.subDict(phases[pi]);

        const dictionary& defaultDict =
            defaults().subDict("materialProperties");

        if (defaultDict.found(phases[pi]+"Material"))
        {
            const dictionary& defaultMatDict =
                defaultDict.subDict(phases[pi]+"Material");
            if (defaultMatDict.found(modelName))
            {
                dict.add(modelName, defaultMatDict.lookup(modelName));
                dict.add(word(modelName+"Coeffs"), dictionary());
                dict.subDict(modelName+"Coeffs").merge
                (
                    defaultMatDict.subOrEmptyDict(modelName+"Coeffs"),
                    // Don't overwrite
                    false
                );
            }
        }
    }
}


void multiphaseMulticomponentState::mergeBinaryPhaseModels
(
    const word& modelName,
    const bool ordered,
    const bool compulsory,
    const word& orderedKeyword
)
{
    dictionary& matDict = system().subDict("materialProperties");
    wordList phases(phaseNames());

    forAll(phases, pi)
    {
        for (label pj = 0; pj < phases.size(); pj++)
        {
            if (pj == pi)
            {
                if (ordered)
                {
                    continue;
                }
                else
                {
                    break;
                }
            }
            // Look for models in the input and the defaults
            // Is there an entry in input?
            phasePairKey key
            (
                phases[pi], phases[pj], ordered, orderedKeyword
            );

            // Check if this model is either unspecified in the input, or
            // is specified in an exact phase key (not a regex). We can then
            // modify the input, but the alternative, where it is specified in
            // a wildcard subdict, cannot be modified since that might influence
            // other phases.

            const dictionary* dict =
                multiphaseThermo::phasePairDictContaining
                (
                    matDict,
                    modelName,
                    phases[pi],
                    phases[pj],
                    ordered
                );
            dictionary* exactDict =
                multiphaseThermo::phasePairDictContaining
                (
                    matDict,
                    modelName,
                    phases[pi],
                    phases[pj],
                    ordered,
                    false
                );

            if (!dict || exactDict)
            {
                // See if a default exists
                const dictionary& defaultBPDict =
                    defaults().subDict("materialProperties");
                const dictionary* defaultDict =
                    multiphaseThermo::phasePairDictContaining
                    (
                        const_cast<dictionary&>(defaultBPDict),
                        modelName,
                        phases[pi],
                        phases[pj],
                        ordered
                    );
                if (defaultDict)
                {
                    // A default does exist. Find or create the dict to add it
                    // to.

                    // Don't use a wildcard entry as that will affect other
                    // phases
                    dictionary* newDict;
                    if (exactDict)
                    {
                        newDict = exactDict;
                    }
                    else
                    {
                        newDict =
                            multiphaseThermo::phasePairDictContaining
                            (
                                matDict,
                                word::null,
                                phases[pi],
                                phases[pj],
                                ordered,
                                false
                            );
                    }
                    if (!newDict)
                    {
                        // Create a new one
                        matDict.add(key.name(), dictionary());
                        newDict = &matDict.subDict(key.name());
                    }
                    newDict->add(modelName, defaultDict->lookup(modelName));
                    newDict->add(word(modelName+"Coeffs"), dictionary());
                    newDict->subDict(modelName+"Coeffs").merge
                    (
                        defaultDict->subOrEmptyDict(modelName+"Coeffs"),
                        // Don't overwrite existing entries
                        false
                    );
                }
                else if (!exactDict && compulsory)
                {
                    FatalIOErrorInFunction(matDict)
                        << modelName << " was not specified"
                        << " in input nor in defaults"
                        << " for phase pair " << key
                        << nl << exit(FatalIOError);
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

multiphaseMulticomponentState::multiphaseMulticomponentState
(
    word region,
    const dictionary& input,
    const dictionary& defaults,
    const stateFunction& master,
    const stateIndex& index,
    word meshName,
    const bool phasicTurbulence
)
:
    turbulenceModelState
    (
        region,
        input,
        defaults,
        master,
        index,
        meshName,
        phasicTurbulence
    ),
    materialNames_(input.lookup("materials"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void multiphaseMulticomponentState::initialise()
{
    // Check that is a material state
    if
    (
        !defaults().lookup<wordList>("materialLibStates").found(stateType())
    )
    {
        FatalErrorInFunction
            << "The state " << stateType() << " is not compatible with the new "
            << "materials library. Valid states are: "
            << defaults().lookup<wordList>("materialLibStates")
            << nl << exit(FatalError);
    }

    // Merge all species data into materialProperties
    initialiseMaterials();

    wordList phases = phaseNames();
    const dictionary& materialProperties =
        system().subDict("materialProperties");
    forAll(phases, phasei)
    {
        const dictionary& phaseDict =
            (
                phases[phasei] == word::null
              ? materialProperties
              : materialProperties.subDict(phases[phasei])
            );
        wordList componentNames =
            phaseDict.lookupOrDefault
            (
                "species", wordList(1, word::null)
            );
        if (componentNames.size() > 1)
        {
            word phaseSuffix;
            if (phases[phasei] != word::null)
            {
                phaseSuffix = word(".") + phases[phasei];
            }

            // set field maps for species
            forAll(componentNames, si)
            {
                // fieldMaps
                stateDict_.subDict("fieldMaps").add
                (
                    word(componentNames[si] + phaseSuffix),
                    input().subDict("fieldMaps").lookup<word>("Ydefault")
                );
            }
            if (phaseSuffix != word::null)
            {
                stateDict_.subDict("fieldMaps").add
                (
                    word("Ydefault" + phaseSuffix),
                    input().subDict("fieldMaps").lookup<word>("Ydefault")
                );
                // Remove old one from input to prevent it being merged later
                input().subDict("fieldMaps").remove("Ydefault");

                // Make the species concentration solvers phasic
                dictionary& fvOpt =
                    input().subDict("system").subDict("fvOptions");
                if (fvOpt.found("speciesConcentrationSolver"))
                {
                    dictionary species =
                        fvOpt.subDict("speciesConcentrationSolver");
                    species.add("phaseName", phases[phasei]);
                    fvOpt.add
                    (
                        word("speciesConcentrationSolver" + phaseSuffix),
                        species
                    );
                    fvOpt.remove("speciesConcentrationSolver");
                }
            }
        }
    }
    turbulenceModelState::initialise();
}

void multiphaseMulticomponentState::finalise()
{
    stateFunction::finalise();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction
} // End namespace Foam

// ************************************************************************* //
