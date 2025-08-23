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
    Foam::saturated

Description
    Model which uses a saturation pressure model for a single species to
    calculate the interface composition.

SourceFiles
    saturated.C

\*---------------------------------------------------------------------------*/

#ifndef saturated_H
#define saturated_H

#include "../eulerianInterfaceCompositionModel/eulerianInterfaceCompositionModel.H"
#include "../../eulerianSaturationModels/eulerianSaturationModel/eulerianSaturationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class eulerianPhasePair;

namespace eulerianInterfaceCompositionModels
{

/*---------------------------------------------------------------------------*\
                         Class saturated Declaration
\*---------------------------------------------------------------------------*/

class saturated
:
    public eulerianInterfaceCompositionModel
{
protected:

    // Private data

        //- saturated species name
        word saturatedName_;

        //- saturated species index
        label saturatedIndex_;

        //- Saturation pressure model
        autoPtr<eulerianSaturationModel> saturationModel_;


    // Private Member Functions

        //- Constant of propotionality between partial pressure and mass
        //  fraction
        tmp<volScalarField> wRatioByP() const;


public:

    //- Runtime type information
    TypeName("saturated");

    // Constructors

        //- Construct from components
        saturated
        (
            const dictionary& dict,
            const eulerianPhasePair& pair
        );


    //- Destructor
    virtual ~saturated();


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
