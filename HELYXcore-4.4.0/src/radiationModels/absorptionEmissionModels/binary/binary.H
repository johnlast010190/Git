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
    (c) 2011-2019 OpenFOAM Foundation

Class
    Foam::radiationModels::absorptionEmissionModels::binary

Description
    Radiation coefficient based on two absorption models

SourceFiles
    binary.C

\*---------------------------------------------------------------------------*/

#ifndef radiationbinary_H
#define radiationbinary_H

#include "absorptionEmissionModels/absorptionEmissionModel/absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{

/*---------------------------------------------------------------------------*\
                           Class binary Declaration
\*---------------------------------------------------------------------------*/

class binary
:
    public absorptionEmissionModel
{
    // Private data

        //- Coefficients dictionary
        dictionary coeffsDict_;

        //- First absorption model
        autoPtr<absorptionEmissionModel> model1_;

        //- Second absorption model
        autoPtr<absorptionEmissionModel> model2_;


public:

    //- Runtime type information
    TypeName("binary");


    // Constructors

        //- Construct from components
        binary
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& typeNameDerived = typeName
        );


    //- Destructor
    virtual ~binary();


    // Member Functions

        // Absorption coefficient

            //- Absorption coefficient for continuous phase
            virtual tmp<volScalarField> aCont(const label bandI = 0) const;

            //- Absorption coefficient for dispersed phase
            virtual tmp<volScalarField> aDisp(const label bandI = 0) const;


        // Emission coefficient

            //- Emission coefficient for continuous phase
            virtual tmp<volScalarField> eCont(const label bandI = 0) const;

            //- Emission coefficient for dispersed phase
            virtual tmp<volScalarField> eDisp(const label bandI = 0) const;


        // Emission contribution

            //- Emission contribution for continuous phase
            virtual tmp<volScalarField> ECont(const label bandI = 0) const;

            //- Emission contribution for continuous phase
            virtual tmp<volScalarField> EDisp(const label bandI = 0) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace absorptionEmissionModels
} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
