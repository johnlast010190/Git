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
    (c) 2011-2020 OpenFOAM Foundation
    (c) 2015 OpenCFD Ltd.

Class
    Foam::radiationModels::absorptionEmissionModels::greyMeanSolid

Description
    greyMeanSolid radiation absorption and emission
    coefficients for solid mixture

SourceFiles
    greyMeanSolid.C

\*---------------------------------------------------------------------------*/

#ifndef greyMeanSolid_H
#define greyMeanSolid_H

#include "absorptionEmissionModels/absorptionEmissionModel/absorptionEmissionModel.H"
#include "solidThermo/solidThermo.H"
#include "mixtures/basicSpecieMixture/basicSpecieMixture.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{

/*---------------------------------------------------------------------------*\
                       Class greyMeanSolid Declaration
\*---------------------------------------------------------------------------*/

class greyMeanSolid
:
    public absorptionEmissionModel
{
private:

    // Private data

        //- Enumering of radiative properties
        enum radiativeProperties
        {
            absorptivity,
            emissivity
        };

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- SLG thermo package
        const solidThermo& thermo_;

        //- Hash table of species names
        HashTable<label> speciesNames_;

        //- Basic multicomponent mixture
        const basicSpecieMixture& mixture_;

        //- List of solid species data
        List<FixedList<scalar, 2>> solidData_;


    // Private member functions

       //- Calculate the volumetric fraction of Yj
       tmp<scalarField> X(const word specie) const;

       //- Calculate the property mixing
       tmp<volScalarField> calc(const label) const;


public:

    //- Runtime type information
    TypeName("greyMeanSolid");


    // Constructors

        //- Construct from components
        greyMeanSolid
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& defaultTypeName = typeName
        );


    //- Destructor
    virtual ~greyMeanSolid();


    // Member Functions

        //- Absorption coefficient for continuous phase
        tmp<volScalarField> aCont(const label bandI = 0) const;

        //- Emission coefficient for continuous phase
        tmp<volScalarField> eCont(const label bandI = 0) const;

        //- Is grey model
        inline bool isGrey() const
        {
            return true;
        }

        //- Helper function to pick up thermo
        const solidThermo& thermo(const fvMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace absorptionEmissionModels
} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
