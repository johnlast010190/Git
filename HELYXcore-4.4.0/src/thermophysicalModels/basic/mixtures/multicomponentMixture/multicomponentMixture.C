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
    (c) 2011-2021 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "mixtures/multicomponentMixture/multicomponentMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
Foam::PtrList<ThermoType>
Foam::multicomponentMixture<ThermoType>::readSpeciesData
(
    const objectRegistry& obr,
    const dictionary& thermoDict
) const
{
    PtrList<ThermoType> speciesData(species_.size());

    forAll(species_, i)
    {
        speciesData.set
        (
            i,
            new ThermoType(obr, thermoDict.subDict(species_[i]))
        );
    }

    return speciesData;
}


template<class ThermoType>
Foam::List<Foam::List<Foam::specieElement>>
Foam::multicomponentMixture<ThermoType>::readSpeciesComposition
(
    const dictionary& thermoDict
) const
{
    List<List<specieElement>> specieCompositions(species_.size());

    // Loop through all species in thermoDict to retrieve
    // the species composition
    forAll(species_, i)
    {
        if (thermoDict.subDict(species_[i]).isDict("elements"))
        {
            const dictionary& elements =
                thermoDict.subDict(species_[i]).subDict("elements");

            const wordList elementsNames(elements.toc());

            specieCompositions[i].resize(elementsNames.size());

            forAll(elementsNames, eni)
            {
                specieCompositions[i][eni].name() = elementsNames[eni];
                specieCompositions[i][eni].nAtoms() =
                    elements.lookupOrDefault(elementsNames[eni], 0);
            }
        }
    }

    return specieCompositions;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::multicomponentMixture<ThermoType>::multicomponentMixture
(
    const dictionary& thermoDict,
    const objectRegistry& obr,
    const word& phaseName
)
:
    basicSpecieMixture
    (
        thermoDict,
        thermoDict.lookup("species"),
        obr,
        phaseName
    ),
    specieThermos_(readSpeciesData(obr, thermoDict)),
    specieCompositions_(readSpeciesComposition(thermoDict))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::multicomponentMixture<ThermoType>::read
(
    const objectRegistry& obr,
    const dictionary& thermoDict
)
{
    specieThermos_ = readSpeciesData(obr, thermoDict);
    specieCompositions_ = readSpeciesComposition(thermoDict);
}


template<class ThermoType>
const Foam::List<Foam::specieElement>&
Foam::multicomponentMixture<ThermoType>::specieComposition
(
    const label speciei
) const
{
    if (specieCompositions_[speciei].empty())
    {
        // Spit an error associated with the lookup of this specie's elements
        refCast<const dictionary>(*this)
            .subDict(species_[speciei])
            .subDict("elements");
    }

    return specieCompositions_[speciei];
}


// ************************************************************************* //
