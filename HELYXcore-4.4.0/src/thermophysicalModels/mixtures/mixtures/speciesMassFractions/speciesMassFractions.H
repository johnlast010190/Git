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
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2024-2025 Engys Ltd.

Class
    Foam::speciesMassFractions

Group
    grpMulticomponentThermophysicalMixtures

Description
    Multi-component mixture.

    Provides a list of mass fraction fields and helper functions to
    query mixture composition.

SourceFiles
    speciesMassFractions.C

\*---------------------------------------------------------------------------*/

#ifndef speciesMassFractions_H
#define speciesMassFractions_H

#include "fields/volFields/volFields.H"
#include "containers/Lists/PtrList/PtrList.H"
#include "db/typeInfo/typeInfo.H"
#include "speciesTable/speciesTable.H"
#include "mixtures/volumeMassFractions/volumeMassFractions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class speciesMassFractions Declaration
\*---------------------------------------------------------------------------*/

class speciesMassFractions
:
    public volumeMassFractions
{
protected:

    // Protected data

        //- Table of specie names
        speciesTable species_;

        //- The name of the default specie
        //  The mass fraction of which is derived from the other species
        //  rather than solved.
        //  Also used as the carrier specie in multicomponent diffusion
        word defaultSpecie_;

        //- The index of the default specie
        label defaultSpecieIndex_;

        //- List of specie active flags
        mutable boolList active_;

        //- Species mass fractions
        PtrList<volScalarField> Y_;


    // Protected Member Functions

        //- Scale the mass fractions to sum to 1
        void correctMassFractions();


public:

    //- Run time type information
    TypeName("speciesMassFractions");


    //- The base class of the mixture
    typedef speciesMassFractions basicMixtureType;


    // Constructors

        //- Read fraction fields
        void createFields
        (
            const dictionary& thermoDict,
            const wordList& specieNames,
            const objectRegistry& obr,
            const word& phaseName
        );

        //- Construct from dictionary, species names, mesh and phase name
        speciesMassFractions
        (
            const dictionary&,
            const wordList& specieNames,
            const objectRegistry&,
            const word&
        );

        //- Construct from dict, object registry and phase name
        speciesMassFractions
        (
            const dictionary& thermoDict,
            const objectRegistry& obr,
            const word& phaseName
        );


    //- Destructor
    virtual ~speciesMassFractions()
    {}


    // Member functions

        //- Make automatic selection of the default specie
        void setDefaultSpecie(const dictionary& dict);

        //- Return the table of species
        inline const speciesTable& species() const;

        //- Does the mixture include this specie?
        inline bool contains(const word& specieName) const;

        //- Return the index of the default specie
        inline label defaultSpecie() const;

        //- Return true for active species
        inline bool active(label speciei) const override;

        //- Return the bool list of active species
        inline const List<bool>& active() const override;

        //- Set speciei active
        inline void setActive(label speciei) const;

        //- Set speciei inactive
        inline void setInactive(label speciei) const;

        //- Return true if the specie should be solved for
        //  i.e. active and not the default specie
        inline bool solve(label speciei) const;

        //- Return the mass-fraction fields
        inline PtrList<volScalarField>& Y();

        //- Return the const mass-fraction fields
        inline const PtrList<volScalarField>& Y() const;

        //- Return the mass/volume fraction fields
        inline PtrList<volScalarField>& fractions() override;

         //- Return the const mass/volume fraction fields
        inline const PtrList<volScalarField>& fractions() const override;

        //- Return the mass-fraction field for a specie given by index
        inline volScalarField& Y(const label i);

        //- Return the const mass-fraction field for a specie given by index
        inline const volScalarField& Y(const label i) const;

        //- Return the mass-fraction field for a specie given by name
        inline volScalarField& Y(const word& specieName);

        //- Return the const mass-fraction field for a specie given by name
        inline const volScalarField& Y(const word& specieName) const;

        //- Return the specie index of the given mass-fraction field
        inline label index(const volScalarField& Yi) const;

        //- Normalise the mass fractions
        //  by clipping positive
        //  and deriving the default specie mass fraction from the other species
        void normalise();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "speciesMassFractionsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
