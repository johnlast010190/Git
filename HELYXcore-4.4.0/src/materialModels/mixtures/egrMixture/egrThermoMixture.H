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
    Foam::egrThermoMixture

Description
    Egr mixture model.

SourceFiles
    egrThermoMixture.C

\*---------------------------------------------------------------------------*/

#ifndef egrThermoMixture_H
#define egrThermoMixture_H

#include "memory/autoPtr/autoPtr.H"
#include "materialModels/materialModel/materialModel.H"
#include "mixtures/speciesMassFractions/speciesMassFractions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class egrThermoMixture Declaration
\*---------------------------------------------------------------------------*/

class egrThermoMixture
:
    public materialModel,
    public egrReactantMixtureModel,
    public egrProductMixtureModel,
    public egrThermoMixtureModel
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

        //- Regress variable
        const volScalarField& b_;

        //- Residual gases
        const volScalarField& egr_;

        //- Stoic ratio
        scalar stoicRatio_;

        //- Harmonic mixing
        bool harmonic_;


public:


    // Constructors

        //- Construct from dictionary
        egrThermoMixture
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& name,
            const word& phaseName,
            const word& specieName
        );

        //- New constructor
        static autoPtr<egrThermoMixture> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );



    //- Destructor
    virtual ~egrThermoMixture()
    {}

    // Member functions

        //- Return the instantiated type name
        TypeName("egrThermoMixture");


    // Member functions

        scalar combustionMixture
        (
            scalar ft,
            scalar b,
            scalar egr,
            scalar oxidant,
            scalar fuel,
            scalar products
        ) const
        {
            if (ft < 0.0001)
            {
                return oxidant;
            }

            scalar fu = b*ft + (1.0 - b)*max(ft - (1.0 - ft)/stoicRatio_, 0.0);
            scalar ox = 1 - ft - (ft - fu)*stoicRatio_;

            fu *= (1.0 - egr);
            ox *= (1.0 - egr);
            if (harmonic_)
            {
                return 1/(fu/fuel + ox/oxidant + (1 - fu - ox)/products);
            }
            return fu*fuel + ox*oxidant + (1 - fu - ox)*products;
        }

        //- The mixture functions are not constant
        bool isegrReactantMixtureConst() const override
        {
            return false;
        }

        //- The mixture functions are not constant
        bool isegrProductMixtureConst() const override
        {
            return false;
        }

        //- The mixture functions are not constant
        bool isegrThermoMixtureConst() const override
        {
            return false;
        }

        declareScalarFieldFunctions(egrReactantMixture);
        declareScalarFieldFunctions(egrProductMixture);
        declareScalarFieldFunctions(egrThermoMixture);

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
