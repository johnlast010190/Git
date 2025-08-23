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
    (c) 2021-2025 Engys Ltd.

Class
    energyInverseTemperature

SourceFiles
    energyInverseTemperature.C

\*---------------------------------------------------------------------------*/

#ifndef energyInverseTemperature_H
#define energyInverseTemperature_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/materialTables/materialTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class energyInverseTemperature Declaration
\*---------------------------------------------------------------------------*/

class energyInverseTemperature
:
    public materialModel,
    public TModel,
    public ThaModel,
    public ThsModel,
    public TeaModel,
    public TesModel
{

    // Static data members

        //- Convergence tolerance of energy -> temperature inversion functions
        static const scalar tol_;

        //- Max number of iterations in energy->temperature inversion functions
        static const int maxIter_;


    // Private data members

        //- Temperature
        volScalarField* TnonConst_;

        //- Reference offset
        scalar Tref_;

        //- Energy
        const volScalarField* he_;

        //- Quick model access
        enum {he, Cpv, limit, ha, hs, ea, es, Cp, Cv, modelsEnumSize_};


    // Private member functions

        //- Lookup energy field name
        word lookupEnergyFieldName();


public:


    // Constructors

        //- Construct from components
        energyInverseTemperature
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    //- Destructor
    virtual ~energyInverseTemperature();


    // Member functions

        //- Return the instantiated type name
        TypeName("energyInverse");

        //- Mixture temperature corresponding to the value of the
        //  thermodynamic property f, given the function f = F(p, T)
        //  and dF(p, T)/dT
        // For cells
        scalar Cell
        (
            const label celli,
            const label heFun,
            const label CpvFun
        ) const;

        //- For patches
        tmp<scalarField> Patch
        (
            const label patchi,
            const label heFun,
            const label CpvFun
        ) const;

        //- For internal fields
        tmp<scalarField> Internal
        (
            const label heFun,
            const label CpvFun
        ) const;

        bool isTConst() const override
        {
            return T_->isConst();
        }

        bool isThsConst() const override
        {
            return T_->isConst();
        }

        bool isThaConst() const override
        {
            return T_->isConst();
        }

        bool isTesConst() const override
        {
            return T_->isConst();
        }

        bool isTeaConst() const override
        {
            return T_->isConst();
        }

        //- Inversion based on selected energy model
        declareScalarFieldFunctions(T);

        //- Inversion with ha
        declareScalarFieldFunctions(Tha);

        //- Inversion with hs
        declareScalarFieldFunctions(Ths);

        //- Inversion with ea
        declareScalarFieldFunctions(Tea);

        //- Inversion with es
        declareScalarFieldFunctions(Tes);

        //- Return model dependence
        materialDepList();

        //- Update table
        virtual void updateTable(const word& modelName) override;

        //- Models casting
        virtual baseModels<scalar>* castScalarModel
        (
            const word& modelName
        ) override;

        //- Update fields
        virtual bool updateScalarField
        (
            const word& fieldName,
            const volScalarField& volField
        ) override;

        //- Re-read data
        virtual bool read() override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
