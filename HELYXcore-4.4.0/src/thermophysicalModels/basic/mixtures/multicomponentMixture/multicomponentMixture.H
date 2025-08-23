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

Class
    Foam::multicomponentMixture

Group
    grpMulticomponentThermophysicalMixtures

Description
    Foam::multicomponentMixture

SourceFiles
    multicomponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentMixture_H
#define multicomponentMixture_H

#include "mixtures/basicSpecieMixture/basicSpecieMixture.H"
#include "containers/HashTables/HashPtrTable/HashPtrTable.H"
#include "specieElement/specieElement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class multicomponentMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class multicomponentMixture
:
    public basicSpecieMixture
{
public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;


private:

    // Private Data

        //- Species data
        PtrList<ThermoType> specieThermos_;

        //- Table of species composition
        List<List<specieElement>> specieCompositions_;


    // Private Member Functions

        //- Read the species data from the given dictionary and return
        PtrList<ThermoType> readSpeciesData
        (
            const objectRegistry& obr,
            const dictionary& thermoDict
        ) const;


        //- Read the species composition from the given dictionary and return
        List<List<specieElement>> readSpeciesComposition
        (
            const dictionary& thermoDict
        ) const;


public:

    // Constructors

        //- Construct from dictionary, obr and phase name
        multicomponentMixture
        (
            const dictionary&,
            const objectRegistry&,
            const word&
        );

        //- Disallow default bitwise copy construction
        multicomponentMixture
        (
            const multicomponentMixture<ThermoType>&
        ) = delete;


    //- Destructor
    virtual ~multicomponentMixture()
    {}


    // Member functions

        //- Return the raw specie thermodynamic data
        inline const PtrList<ThermoType>& specieThermos() const
        {
            return specieThermos_;
        }

        //- Read dictionary
        void read(const objectRegistry&, const dictionary&);

        //- Return thermo based on index
        inline const ThermoType& specieThermo(const label speciei) const
        {
            return specieThermos_[speciei];
        }

        //- Return composition based on index
        const List<specieElement>& specieComposition
        (
            const label speciei
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mixtures/multicomponentMixture/multicomponentMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
