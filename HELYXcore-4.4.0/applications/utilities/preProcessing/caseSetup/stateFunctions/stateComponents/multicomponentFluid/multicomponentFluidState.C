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

#include "multicomponentFluidState.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "stateComponents/multiphaseMulticomponent/multiphaseMulticomponentState.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace stateFunctions
{

defineTypeNameAndDebug(multicomponentFluidState, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

dictionary multicomponentFluidState::mixtureProperties
(
    const word& regionName,
    const word& inertName,
    const dictionary& mat
)
{
    dictionary matDict;

    matDict.add("inertSpecie", inertName);

    // select thermoType based on inert
    matDict.add
    (
        "thermoType",
        mat.subDict(inertName).subDict("thermoType")
    );

    // default chemical species thermo data
    matDict.add("#include", fileName("speciesProperties"));

    return matDict;
}


dictionary multicomponentFluidState::thermoSpeciesProperties
(
    const wordList& componentNames,
    const dictionary& mat
)
{
    dictionary specDict;
    forAll(componentNames, specI)
    {
        specDict.add(componentNames[specI], dictionary());
        dictionary& species(specDict.subDict(componentNames[specI]));
        const dictionary& matSpecies = mat.subDict(componentNames[specI]);
        //add mixture properties (for each species)
        species.add("specie", matSpecies.subDict("specie"));
        if (matSpecies.found("elements"))
        {
            species.add("elements", matSpecies.subDict("elements"));
        }
        species.add("thermodynamics", matSpecies.subDict("thermodynamics"));
        species.add("transport", matSpecies.subDict("transport"));
        species.add
        (
            "equationOfState",
            matSpecies.subOrEmptyDict("equationOfState")
        );
    }

    // TODO make sure all the thermo models between species are the
    // same and we have enough info
    return specDict;
}


dictionary multicomponentFluidState::materialSpeciesProperties
(
    const wordList& componentNames,
    const dictionary& mat
)
{
    dictionary specDict;
    forAll(componentNames, specI)
    {
        specDict.add(componentNames[specI], dictionary());
        dictionary& species(specDict.subDict(componentNames[specI]));
        const dictionary& matSpecies = mat.subDict(componentNames[specI]);
        species.merge
        (
            multiphaseMulticomponentState::speciesProperties
            (
                mat.lookup<word>("materialType"), matSpecies
            )
        );
    }

    return specDict;
}


dictionary multicomponentFluidState::combustionProperties
(
    const wordList& componentNames,
    const dictionary& combustionInput
)
{
    dictionary combustionDict(combustionInput);
    combustionDict.add("combustionModel", "none");
    return combustionDict;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

multicomponentFluidState::multicomponentFluidState
(
    word region,
    const dictionary& input,
    const dictionary& defaults,
    const stateFunction& master,
    const stateIndex& index,
    word meshName
)
:
    singlePhaseFluidState
    (
        region,
        input,
        defaults,
        master,
        index,
        meshName
    ),
    region_(region)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void multicomponentFluidState::initialiseMaterials(dictionary& materialProperties)
{
    // Merge all species thermo data into materialProperties
    const dictionary& inDict = input().subDict("materialProperties");
    materialProperties.add("materialType", inDict.lookup("materialType"));
    materialProperties.add("mixture", inDict.lookupOrDefault<word>("mixture", "standardMixture"));
    materialProperties.add("energy", inDict.lookupOrDefault<word>("energy", "sensibleEnthalpy"));

    if (input().lookup<wordList>("state").found("buoyant") || inDict.found("buoyancyModel"))
    {
        if (inDict.found("buoyancyModel"))
        {
            materialProperties.add
            (
                "buoyancyModel",
                inDict.lookup<word>("buoyancyModel")
            );

            if (inDict.found("buoyancyModelCoeffs"))
            {
                materialProperties.add
                (
                    "buoyancyModelCoeffs",
                    inDict.subDict("buoyancyModelCoeffs")
                );
            }
        }
        else if (input().isDict("buoyantSteady"))
        {
            materialProperties.add
            (
                "buoyancyModel",
                input().subDict("buoyantSteady").subDict("materialProperties").lookupOrDefault<word>
                (
                    "buoyancyModel",
                    "rhoModel"
                )
            );
        }
        else
        {
            materialProperties.add("buoyancyModel", "rhoModel");
        }
    }

    materialProperties.merge(defaults().subDict("materialProperties").subDict("referenceValues"));
    if (inDict.found("referenceFields"))
    {
        materialProperties.subDict("referenceFields").merge
        (
            inDict.subOrEmptyDict("referenceFields")
        );
    }

    forAll(componentNames(), specI)
    {
        materialProperties.add(componentNames()[specI], dictionary(), false);
        materialProperties.subDict(componentNames()[specI]).merge
        (
            defaults().subDict("materialProperties").subOrEmptyDict
            (
                componentNames()[specI] + "Material"
            )
        );

        if (input().found("materialProperties"))
        {
            materialProperties.subDict(componentNames()[specI]).merge
            (
                input().subDict("materialProperties").subOrEmptyDict(componentNames()[specI])
            );
        }
    }

    if (compressibility() == ctComp)
    {
        system().add("materialProperties", dictionary(), false);
        system().subDict("materialProperties").add("species", componentNames());

        dictionary mixtureDict;
        if (input().subDict("materialProperties").found("inertSpecie"))
        {
            mixtureDict.add
            (
                "inertSpecie",
                input().subDict("materialProperties").lookup<word>("inertSpecie")
            );
        }
        mixtureDict.add("materialType", materialProperties.lookup<word>("materialType"));
        mixtureDict.add("mixture", materialProperties.lookupOrDefault<word>("mixture", "standardMixture"));
        mixtureDict.add("energy", materialProperties.lookupOrDefault<word>("energy", "sensibleEnthalpy"));
        if (materialProperties.found("buoyancyModel"))
        {
            mixtureDict.add("buoyancyModel", materialProperties.lookup<word>("buoyancyModel"));
            if (materialProperties.found("buoyancyModelCoeffs"))
            {
                mixtureDict.add
                (
                    "buoyancyModelCoeffs",
                    materialProperties.subDict("buoyancyModelCoeffs")
                );
            }
        }
        mixtureDict.add("referenceFields", materialProperties.subOrEmptyDict("referenceFields"));
        // multi-species case
        // use thermo properties of inert (last species in list)
        system().subDict("materialProperties").merge(mixtureDict);

        // species thermo data
        system().subDict("materialProperties").merge
        (
            materialSpeciesProperties(componentNames(), materialProperties)
        );
    }
}


void multicomponentFluidState::initialiseThermo(dictionary& materialProperties)
{
    forAll(componentNames(), specI)
    {
        materialProperties.add(componentNames()[specI], dictionary(), false);
        if (input().found("materialProperties"))
        {
            materialProperties.subDict(componentNames()[specI]).merge
            (
                input().subDict("materialProperties").subOrEmptyDict(componentNames()[specI])
            );
        }
    }

    if (compressibility() == ctComp)
    {
        const dictionary mixtureDict
        (
            mixtureProperties
            (
                region_,
                componentNames().last(),
                materialProperties
            )
        );
        constant().add("thermophysicalProperties", dictionary(), false);

        // multi-species case
        // use thermo properties of inert (last species in list)
        constant().subDict("thermophysicalProperties").merge(mixtureDict);

        // species thermo data
        constant().add("speciesProperties", dictionary(), false);
        constant().subDict("speciesProperties").merge
        (
            thermoSpeciesProperties(componentNames(), materialProperties)
        );

        constant().subDict("speciesProperties").add
        (
            "species",
            componentNames()
        );

        // TODO: Reaction properties have to be added to chemistryProperties
        //       dict in the future
        // constant().add("reactionProperties", dictionary(), false);
        // dictionary reactionInput;
        // reactionInput.merge
        // (
        //     defaults().subOrEmptyDict("constant").subOrEmptyDict
        //     (
        //         "reactionProperties"
        //     )
        // );
        // reactionInput.merge
        // (
        //     input().subOrEmptyDict("constant").subOrEmptyDict
        //     (
        //         "reactionProperties"
        //     )
        // );

        // // Empty reactions for now
        // dictionary reactionDict(reactionInput);
        // reactionDict.add("reactions", dictionary());
        // constant().subDict("reactionProperties").merge(reactionDict);
    }
}


void multicomponentFluidState::initialise()
{
    const bool isMaterials =
        defaults().lookup<wordList>("materialLibStates").found(stateType());
    if (materialName() != word::null)
    {
        singlePhaseFluidState::initialise();
    }
    else
    {
        // Material properties
        dictionary materialProperties;

        // Merge all species thermo data into materialProperties
        if (isMaterials)
        {
            initialiseMaterials(materialProperties);
        }
        else
        {
            initialiseThermo(materialProperties);
        }
        if (compressibility() == ctComp)
        {
            // generic combustion setup
            constant().add("combustionProperties", dictionary(), false);
            dictionary combustionInput;
            combustionInput.merge
            (
                defaults().subOrEmptyDict("constant").subOrEmptyDict
                (
                    "combustionProperties"
                )
            );
            combustionInput.merge
            (
                input().subOrEmptyDict("constant").subOrEmptyDict
                (
                    "combustionProperties"
                )
            );
            constant().subDict("combustionProperties").merge
            (
                multicomponentFluidState::combustionProperties
                (
                    componentNames(), combustionInput
                )
            );
            // set field maps for species
            forAll(componentNames(), si)
            {
                // fieldMaps
                stateDict_.subDict("fieldMaps").add
                (
                    word(componentNames()[si]),
                    input().subDict("fieldMaps").lookup("Ydefault")
                );
            }
        }
        else
        {
            FatalErrorInFunction
                << "Invalid compressibility type: "
                << compressibility() <<" for this state."
                << exit(FatalError);
        }
    }
    turbulenceModelState::initialise();
}


void multicomponentFluidState::finalise()
{
    stateFunction::finalise();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction
} // End namespace Foam

// ************************************************************************* //
