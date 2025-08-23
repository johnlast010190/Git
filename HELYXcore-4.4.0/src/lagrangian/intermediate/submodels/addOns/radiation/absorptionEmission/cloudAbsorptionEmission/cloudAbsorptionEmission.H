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
    Foam::radiationModels:::cloudAbsorptionEmission

Description
    Retrieves absorption/emission data from a cloud object

SourceFiles
    cloudAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef radiationCloudAbsorptionEmission_H
#define radiationCloudAbsorptionEmission_H

#include "absorptionEmissionModels/absorptionEmissionModel/absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{

/*---------------------------------------------------------------------------*\
                  Class cloudAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class cloudAbsorptionEmission
:
    public absorptionEmissionModel
{
    // Private data

        //- Coefficients dictionary
        dictionary coeffsDict_;

        //- Cloud name(s)
        const wordList cloudNames_;


public:

    //- Runtime type information
    TypeName("cloud");


    // Constructors

        //- Construct from components
        cloudAbsorptionEmission
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& derivedClassName = typeName
        );


    //- Destructor
    virtual ~cloudAbsorptionEmission();


    // Member Operators

        //- Absorption coefficient for dispersed phase
        virtual tmp<volScalarField> aDisp(const label bandI = 0) const;

        //- Emission coefficient for dispersed phase
        virtual tmp<volScalarField> eDisp(const label bandI = 0) const;

        //- Return emission contribution for dispersed phase
        virtual tmp<volScalarField> EDisp(const label bandI = 0) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace absorptionEmissionModels
} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
