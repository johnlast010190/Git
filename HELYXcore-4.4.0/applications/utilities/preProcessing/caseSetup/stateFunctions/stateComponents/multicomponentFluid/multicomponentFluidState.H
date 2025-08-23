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

Class
    Foam::multicomponentFluidState

Description
    Abstract, multicomponentFluid state transport/thermo assembly

    This has been replaced by multiphaseMulticomponentFluidState for USF states
    and should be phased out as legacy solvers are removed.

SourceFiles
    multicomponentFluidState.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentFluidState_H
#define multicomponentFluidState_H

#include "stateComponents/singlePhaseFluid/singlePhaseFluidState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class multicomponentFluidState Declaration
\*---------------------------------------------------------------------------*/

class multicomponentFluidState
:
    public singlePhaseFluidState
{

protected:

    // Protected Data
    word region_;


public:

    // Static member functions

        //- Single phase multi-species compressible thermophyscial properties
        static dictionary mixtureProperties
        (
            const word& regionName,
            const word& inertSpecie,
            const dictionary& mat
        );

        //- Single phase multi-species compressible species properties
        static dictionary thermoSpeciesProperties
        (
            const wordList& specieNames,
            const dictionary& mat
        );

        //- Single phase multi-species compressible material properties
        static dictionary materialSpeciesProperties
        (
            const wordList& specieNames,
            const dictionary& mat
        );

        //- Multi-species combustion properties
        static dictionary combustionProperties
        (
            const wordList& speciesNames,
            const dictionary& combustionInput
        );

    //- Runtime type information
    TypeName("multicomponentFluid");

    // Constructors

        //- Construct from components
        multicomponentFluidState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );


    //- Destructor
    virtual ~multicomponentFluidState()
    {}


    // Member Functions

        //- Assemble field definitions
        virtual void initialise();

        //- Initialise material library
        void initialiseMaterials(dictionary& materialProperties);

        //- Initialise thermo/transport library
        void initialiseThermo(dictionary& materialProperties);

        //- Final manipulations
        virtual void finalise();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
