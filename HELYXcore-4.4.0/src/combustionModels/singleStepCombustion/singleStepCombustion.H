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
    (c) 2024 Engys Ltd.

Class
    Foam::combustionModels::singleStepCombustion

Group
    grpCombustionModels

Description
    Base class for combustion models using basicSpecieMixture.

SourceFiles
    singleStepCombustion.C

\*---------------------------------------------------------------------------*/

#ifndef singleStepCombustion_H
#define singleStepCombustion_H

#include "mixtures/basicSpecieMixture/basicSpecieMixture.H"
#include "reaction/reaction/reaction.H"
#include "combustionModel/combustionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                    Class singleStepCombustion Declaration
\*---------------------------------------------------------------------------*/

class singleStepCombustion
:
    public combustionModel
{
protected:

    // Protected data

        //- Reference to the mixture
        const basicSpecieMixture& mixture_;

        //- The single-step reaction
        reaction reaction_;

        //- Stoichiometric air-fuel mass ratio
        dimensionedScalar stoicRatio_;

        //- Stoichiometric oxygen-fuel mass ratio
        dimensionedScalar s_;

        //- Heat of combustion [J/Kg]
        dimensionedScalar qFuel_;

        //- Stoichiometric coefficient for the reaction.
        scalarList specieStoichCoeffs_;

        //- Mass concentrations at stoichiometric mixture for fres.
        scalarList Yprod0_;

        //- List of components residual
        PtrList<volScalarField> fres_;

        //- Fuel specie index
        label fuelIndex_;

        //- List to indicate if specie is produced/consumed
        List<int> specieProd_;

        //- Fuel consumption rate
        volScalarField wFuel_;

        //- Semi-implicit (true) or explicit (false) treatment
        bool semiImplicit_;


    // Protected member functions

        //- Calculate qFuel
        void calculateqFuel();

        //- Calculate air/fuel and oxygen/fuel ratio
        void massAndAirStoichRatios();

        //- Calculate maximum products at stoichiometric mixture
        void calculateMaxProducts();


public:

    // Constructors

        //- Construct from components
        singleStepCombustion
        (
            const word& modelType,
            const fluidMulticomponentThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
        );

        //- Disallow default bitwise copy construction
        singleStepCombustion(const singleStepCombustion&) = delete;


    //- Destructor
    virtual ~singleStepCombustion();


    // Member Functions

        // Access functions

            //- Return the stoichiometric air-fuel mass ratio
            inline const dimensionedScalar& stoicRatio() const;

            //- Return the Stoichiometric oxygen-fuel mass ratio
            inline const dimensionedScalar& s() const;

            //- Return the heat of combustion [J/Kg]
            inline const dimensionedScalar& qFuel() const;

            //- Return the stoichiometric coefficient for the reaction
            inline const List<scalar>& specieStoichCoeffs() const;

            //- Return the list of components residual
            inline tmp<volScalarField> fres(const label index) const;

            //- Return the fuel specie index
            inline label fuelIndex() const;

            //- Return the list to indicate if specie is produced/consumed
            inline const List<int>& specieProd() const;

            //- Return the list of products mass concentrations
            inline const scalarList& Yprod0() const;

            //- Reactions
            inline const reaction& singleReaction() const;

            //- Return the mixture
            inline const basicSpecieMixture& mixture() const;


        //- Fuel consumption rate matrix
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

        //- Heat release rate [kg/m/s3]
        virtual tmp<volScalarField> Qdot() const;

        //- Update properties from given dictionary
        virtual bool read();

        //- Calculates the residual for all components
        void fresCorrect();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const singleStepCombustion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "singleStepCombustion/singleStepCombustionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
