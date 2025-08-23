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
    Foam::multiphaseMulticomponentState

Description
    Base state for new material library which can handle one or more phases
    each with one or more components

SourceFiles
    multiphaseMulticomponentState.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseMulticomponentState_H
#define multiphaseMulticomponentState_H

#include "stateComponents/turbulenceModel/turbulenceModelState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class multiphaseMulticomponentState Declaration
\*---------------------------------------------------------------------------*/

class multiphaseMulticomponentState
:
    public turbulenceModelState
{
    // Private Data

        //- The material names as specified in caseSetupDict. Either the phase
        //  materials (for multiphase) or species materials (for single phase)
        //  depending on the derived state
        wordList materialNames_;


protected:

    // Protected member functions

        //- Initialise material properties
        void initialiseMaterials();

        //- Initialise phase sub-dict of material properties
        void initialisePhase
        (
            const word& phaseName,
            dictionary& materialProperties,
            const dictionary& inDict
        );

        //- Initialise component sub-dict of material properties
        void initialiseComponent
        (
            const word& componentName,
            dictionary& materialProperties,
            const dictionary& inDict,
            const dictionary& inDictParent
        );

        //- translate generic phase field gradient entry to phase specific entry
        void addAlphaTypeGrad(const word& phaseName);

        //- Merge defaults for stipulated models in each phase dict
        void mergePhaseModels(const word& modelName);

        //- Merge defaults for stipulated binary models, for each phase pair
        void mergeBinaryPhaseModels
        (
            const word& modelName,
            const bool ordered,
            const bool compulsory=true,
            const word& orderedKeyword="in"
        );

public:

    // Static member functions

        //- Species properties
        static dictionary speciesProperties
        (
            const word& matType,
            const dictionary& matSpecies
        );

        //- Runtime type information
        TypeName("multiphaseMulticomponent");


    // Constructors

        //- Construct from components
        multiphaseMulticomponentState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null,
            const bool phasicTurbulence = false
        );

        //- Disallow default bitwise copy construction
        multiphaseMulticomponentState
        (
            const multiphaseMulticomponentState&
        ) = delete;


    //- Destructor
    virtual ~multiphaseMulticomponentState()
    {}


    // Member Functions

        // Access

            //- The list of materials as specified in caseSetupDict - these
            //  are the top-level materials, either species in a single-phase
            //  state or phase materials in a multi-phase state
            virtual wordList materialNames() const
            {
                return wordList(materialNames_);
            }

            //- Override to provide the names of the phases
            virtual wordList phaseNames() const override
            {
                // Default is the single-phase case unless
                // 'materialType multiphase' is specified in the input
                if
                (
                    input().subOrEmptyDict
                    (
                        "materialProperties"
                    ).lookupOrDefault
                    (
                        "materialType", word::null
                    ) == "multiphase"
                )
                {
                    return materialNames();
                }
                else
                {
                    return wordList(1, word::null);
                }
            }

        //- Assemble field definitions
        virtual void initialise() override;

        //- Final manipulations
        virtual void finalise() override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multiphaseMulticomponentState&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
