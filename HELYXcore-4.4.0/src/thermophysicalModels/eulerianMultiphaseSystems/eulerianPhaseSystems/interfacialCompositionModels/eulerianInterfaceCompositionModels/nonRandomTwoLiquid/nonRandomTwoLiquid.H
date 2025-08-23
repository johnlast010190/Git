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
    Foam::nonRandomTwoLiquid

Description
    Non ideal law for the mixing of two species. A separate composition model
    is given for each species. The composition of a species is equal to the
    value given by the model, scaled by the species fraction in the bulk of the
    other phase, and multiplied by the activity coefficient for that species.
    The gas behaviour is assumed ideal; i.e. the fugacity coefficient is taken
    as equal to 1.

SourceFiles
    nonRandomTwoLiquid.C

\*---------------------------------------------------------------------------*/

#ifndef nonRandomTwoLiquid_H
#define nonRandomTwoLiquid_H

#include "../eulerianInterfaceCompositionModel/eulerianInterfaceCompositionModel.H"
#include "../../eulerianSaturationModels/eulerianSaturationModel/eulerianSaturationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class eulerianPhasePair;

namespace eulerianInterfaceCompositionModels
{

/*---------------------------------------------------------------------------*\
                         Class nonRandomTwoLiquid Declaration
\*---------------------------------------------------------------------------*/

class nonRandomTwoLiquid
:
    public eulerianInterfaceCompositionModel
{
    // Private data

        //- Activity coefficient for species 1
        volScalarField gamma1_;

        //- Activity coefficient for species 2
        volScalarField gamma2_;

        //- Name of species 1
        word species1Name_;

        //- Name of species 2
        word species2Name_;

        //- Indiex of species 1 within this thermo
        label species1Index_;

        //- Indiex of species 2 within this thermo
        label species2Index_;

        //- Non-randomness constant parameter for species 1
        dimensionedScalar alpha12_;

        //- Non-randomness constant parameter for species 2
        dimensionedScalar alpha21_;

        //- Non-randomness linear paramater for species 1
        dimensionedScalar beta12_;

        //- Non-randomness linear paramater for species 2
        dimensionedScalar beta21_;

        //- Interaction parameter model for species 1
        autoPtr<eulerianSaturationModel> saturationModel12_;

        //- Interaction parameter model for species 2
        autoPtr<eulerianSaturationModel> saturationModel21_;

        //- Composition model for species 1
        autoPtr<eulerianInterfaceCompositionModel> speciesModel1_;

        //- Composition model for species 2
        autoPtr<eulerianInterfaceCompositionModel> speciesModel2_;


public:

    //- Runtime type information
    TypeName("nonRandomTwoLiquid");

    // Constructors

        //- Construct from components
        nonRandomTwoLiquid
        (
            const dictionary& dict,
            const eulerianPhasePair& pair
        );


    //- Destructor
    virtual ~nonRandomTwoLiquid();


    // Member Functions

        //- Update the composition
        virtual void update(const volScalarField& Tf);

        //- The interface species fraction
        virtual tmp<volScalarField> Yf
        (
            const word& speciesName,
            const volScalarField& Tf
        ) const;

        //- The interface species fraction derivative w.r.t. temperature
        virtual tmp<volScalarField> YfPrime
        (
            const word& speciesName,
            const volScalarField& Tf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace eulerianInterfaceCompositionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
