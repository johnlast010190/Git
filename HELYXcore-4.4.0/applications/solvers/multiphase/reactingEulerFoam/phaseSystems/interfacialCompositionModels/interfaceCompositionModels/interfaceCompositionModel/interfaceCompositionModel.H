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
    (c) 2015-2019 OpenFOAM Foundation

Class
    Foam::interfaceCompositionModel

Description
    Generic base class for interface composition models. These models describe
    the composition in phase 1 of the supplied pair at the interface with phase
    2.

SourceFiles
    interfaceCompositionModel.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceCompositionModel_H
#define interfaceCompositionModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fields/volFields/volFields.H"
#include "db/dictionary/dictionary.H"
#include "primitives/strings/lists/hashedWordList.H"
#include "rhoMulticomponentThermo/rhoMulticomponentThermo.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

namespace Foam
{

class phaseModel;
class phasePair;

/*---------------------------------------------------------------------------*\
                Class interfaceCompositionModel Declaration
\*---------------------------------------------------------------------------*/

class interfaceCompositionModel
{
protected:

    // Private Data

        //- Phase pair
        const phasePair& pair_;

        //- Names of the transferring species
        const hashedWordList species_;

        //- Lewis number
        const dimensionedScalar Le_;

        //- Multi-component thermo model for this side of the interface
        const rhoMulticomponentThermo& thermo_;

        //- General thermo model for the other side of the interface
        const rhoThermo& otherThermo_;


public:

    //- Runtime type information
    TypeName("interfaceCompositionModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            interfaceCompositionModel,
            dictionary,
            (
                const dictionary& dict,
                const phasePair& pair
            ),
            (dict, pair)
        );


    // Constructors

        //- Construct from a dictionary and a phase pair
        interfaceCompositionModel
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~interfaceCompositionModel();


    // Selectors

        static autoPtr<interfaceCompositionModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions

        // Access

            //- Return the phase pair
            inline const phasePair& pair() const;

            //- Return the transferring species names
            inline const hashedWordList& species() const;

            //- Return the thermo
            inline const rhoMulticomponentThermo& thermo() const;

            //- Return the composition
            inline const basicSpecieMixture& composition() const;

            //- Return the other thermo
            inline const rhoThermo& otherThermo() const;

            //- Return whether the other side has a multi-specie composition
            inline bool otherHasComposition() const;

            //- Return the other composition
            inline const basicSpecieMixture& otherComposition() const;


        // Evaluation

            //- Interface mass fraction
            virtual tmp<volScalarField> Yf
            (
                const word& speciesName,
                const volScalarField& Tf
            ) const = 0;

            //- The interface mass fraction derivative w.r.t. temperature
            virtual tmp<volScalarField> YfPrime
            (
                const word& speciesName,
                const volScalarField& Tf
            ) const = 0;

            //- Mass fraction difference between the interface and the field
            virtual tmp<volScalarField> dY
            (
                const word& speciesName,
                const volScalarField& Tf
            ) const;

            //- Mass diffusivity
            virtual tmp<volScalarField> D
            (
                const word& speciesName
            ) const;

            //- Latent heat
            virtual tmp<volScalarField> L
            (
                const word& speciesName,
                const volScalarField& Tf
            ) const;

            //- Add latent heat flow rate to total
            virtual void addMDotL
            (
                const volScalarField& K,
                const volScalarField& Tf,
                volScalarField& mDotL,
                volScalarField& mDotLPrime
            ) const;


        //- Update the composition
        virtual void update(const volScalarField& Tf) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interfaceCompositionModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
