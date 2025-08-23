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
    (c) 2025 Engys Ltd.

Class
    Foam::veryInhomogeneousThermoMixture

Description
    Very inhomogeneous mixture model.

SourceFiles
    veryInhomogeneousThermoMixture.C

\*---------------------------------------------------------------------------*/

#ifndef veryInhomogeneousThermoMixture_H
#define veryInhomogeneousThermoMixture_H

#include "memory/autoPtr/autoPtr.H"
#include "materialModels/materialModel/materialModel.H"
#include "mixtures/speciesMassFractions/speciesMassFractions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class veryInhomogeneousThermoMixture Declaration
\*---------------------------------------------------------------------------*/

class veryInhomogeneousThermoMixture
:
    public materialModel,
    public veryInhomogeneousReactantMixtureModel,
    public veryInhomogeneousProductMixtureModel,
    public veryInhomogeneousThermoMixtureModel
{
    // Private member data

        //- Quick model access
        enum {fuelInd, oxidantInd, productsInd, modelsEnumSize_};

        //- Reference to the class providing mixing fractions
        const speciesMassFractions& frac_;

        //- Lookup or construct fractions class
        const speciesMassFractions& lookupOrConstructBase();

        //- Mixture fraction
        const volScalarField& ft_;

        //- Fuel mass fraction
        const volScalarField& fu_;

        //- Regress variable
        const volScalarField& b_;

        //- Stoic ratio
        scalar stoicRatio_;

        scalar fres(scalar ft) const
        {
            return max(ft - (1.0 - ft)/stoicRatio_, 0.0);
        }

        //- Harmonic mixing
        bool harmonic_;


public:


    // Constructors

        //- Construct from dictionary
        veryInhomogeneousThermoMixture
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& name,
            const word& phaseName,
            const word& specieName
        );

        //- New constructor
        static autoPtr<veryInhomogeneousThermoMixture> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );



    //- Destructor
    virtual ~veryInhomogeneousThermoMixture()
    {}

    // Member functions

        //- Return the instantiated type name
        TypeName("veryInhomogeneousThermoMixture");


    // Member functions

        scalar combustionMixture
        (
            scalar ft,
            scalar fu,
            scalar oxidant,
            scalar fuel,
            scalar products
        ) const
        {
            if (ft < 0.0001)
            {
                return oxidant;
            }
            else
            {
                scalar ox = 1 - ft - (ft - fu)*stoicRatio_;
                scalar pr = 1 - fu - ox;
                if (harmonic_)
                {
                    return 1.0/(fu/fuel + ox/oxidant + pr/products);
                }
                return fu*fuel + ox*oxidant + pr*products;
            }
        }

        //- The mixture functions are not constant
        bool isveryInhomogeneousReactantMixtureConst() const override
        {
            return false;
        }

        //- The mixture functions are not constant
        bool isveryInhomogeneousProductMixtureConst() const override
        {
            return false;
        }

        //- The mixture functions are not constant
        bool isveryInhomogeneousThermoMixtureConst() const override
        {
            return false;
        }

        declareScalarFieldFunctions(veryInhomogeneousReactantMixture);
        declareScalarFieldFunctions(veryInhomogeneousProductMixture);
        declareScalarFieldFunctions(veryInhomogeneousThermoMixture);

        //- Return model dependence
        materialDepList();

        //- Update table
        virtual void updateTable(const word& modelName) override;

        //- Cast to scalar model
        virtual baseModels<scalar>* castScalarModel
        (
            const word& modelName
        ) override;

        //- Is the model a mixture
        virtual bool mixture() const override
        {
            return true;
        }

        //- Re-read data
        virtual bool read() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
