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
    Foam::radiationModels::absorptionEmission::greyMeanCombustion

Description
    greyMeanCombustion radiation absorption and emission coefficients for
    continuous phase. Exactly as greyMean, but with additional effects of the
    combustion heat release rate included.

See also
    Foam::radiationModels::absorptionEmissionModels::greyMean

SourceFiles
    greyMeanCombustion.C

\*---------------------------------------------------------------------------*/

#ifndef greyMeanCombustion_H
#define greyMeanCombustion_H

#include "interpolations/interpolationLookUpTable/interpolationLookUpTable.H"
#include "absorptionEmissionModels/absorptionEmissionModel/absorptionEmissionModel.H"
#include "absorptionEmissionModels/greyMean/greyMean.H"
#include "containers/HashTables/HashTable/HashTable.H"
#include "radiationModels/fvDOM/absorptionCoeffs/absorptionCoeffs.H"
#include "fluidThermo/fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{

/*---------------------------------------------------------------------------*\
                         Class greyMeanCombustion Declaration
\*---------------------------------------------------------------------------*/

class greyMeanCombustion
:
    public greyMean
{
    // Private data

        //- Emission constant coefficient
        const scalar EhrrCoeff_;


public:

    //- Runtime type information
    TypeName("greyMeanCombustion");


    // Constructors

        //- Construct from components
        greyMeanCombustion
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& typeNameDerived = typeName
        );


    //- Destructor
    virtual ~greyMeanCombustion();


    // Member Functions

        //- Emission contribution for continuous phase
        tmp<volScalarField> ECont(const label bandI = 0) const;

        inline bool isGrey() const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace absorptionEmissionModels
} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
