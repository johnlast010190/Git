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
    Foam::singlePhaseFluidState

Description
    Abstract, singlePhaseFluid state transport/thermo assembly

SourceFiles
    singlePhaseFluidState.C

\*---------------------------------------------------------------------------*/

#ifndef singlePhaseFluidState_H
#define singlePhaseFluidState_H

#include "stateComponents/turbulenceModel/turbulenceModelState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class singlePhaseFluidState Declaration
\*---------------------------------------------------------------------------*/

class singlePhaseFluidState
:
    public turbulenceModelState
{


protected:
    // Protected Data


    // Protected member functions


private:

    // Private Data

        word materialName_;

        wordList componentNames_;



public:

    //- static member functions

        //- single phase incompressible transport properties
        static dictionary incompressibleProperties
        (
            const word& materialName,
            const dictionary& materialProperties
        );

        //- single phase compressible thermophyscial properties
        static dictionary compressibleProperties
        (
            const word& materialName,
            const dictionary& materialProperties
        );


    //- Runtime type information
    TypeName("singlePhaseFluid");

    // Constructors

        //- Construct from components
        singlePhaseFluidState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );

        //- Disallow default bitwise copy construction
        singlePhaseFluidState(const singlePhaseFluidState&) = delete;


    //- Destructor
    virtual ~singlePhaseFluidState()
    {}


    // Member Functions

        // Access

            //const reference to name of fluid
            inline const word& materialName() const
            {
                return materialName_;
            }

            //const reference to name of fluid
            inline const wordList& componentNames() const
            {
                return componentNames_;
            }


        // Edit

            //- assemble field definitions
            virtual void initialise();

            //- single phase compressible thermophyscial properties
            dictionary setupMaterialProperties
            (
                const word& materialName,
                const dictionary& materialProperties
            );

            //- Initialise material library
            void initialiseMaterials
            (
                dictionary& materialProperties,
                const dictionary& matDefaults,
                const dictionary& matInputs,
                const bool isMatInInput
            );

            //- Initialise thermo library
            void initialiseThermo
            (
                dictionary& materialProperties,
                const dictionary& matDefaults,
                const dictionary& matInputs,
                const bool isMatInInput
            );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const singlePhaseFluidState&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
