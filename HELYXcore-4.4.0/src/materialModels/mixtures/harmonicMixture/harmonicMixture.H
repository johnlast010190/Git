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
    (c) 2021-2024 Engys Ltd.

Class
    Foam::harmonicMixture

Description
    Homogeneus mixture model.

SourceFiles
    harmonicMixture.C

\*---------------------------------------------------------------------------*/

#ifndef harmonicMixture_H
#define harmonicMixture_H

#include "memory/autoPtr/autoPtr.H"
#include "materialModels/materialModel/materialModel.H"
#include "mixtures/speciesMassFractions/speciesMassFractions.H"
#include "mixtures/phaseVolumeFractions/phaseVolumeFractions.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class harmonicMixture Declaration
\*---------------------------------------------------------------------------*/

class harmonicMixture
:
    public materialModel,
    public harmonicMixtureModel
{
    // Private member data

        //- Mass fractions class
        const volumeMassFractions& frac_;

        //- Lookup or construct mass fraction class
        const volumeMassFractions& lookupOrConstructBase();

public:

    // Constructors

        //- Construct from dictionary
        harmonicMixture
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& name,
            const word& phaseName,
            const word& specieName
        );

        //- New constructor
        static autoPtr<harmonicMixture> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );



    //- Destructor
    virtual ~harmonicMixture()
    {}

    // Member functions

        //- Return the instantiated type name
        TypeName("harmonicMixture");


    // Member functions

        bool isharmonicMixtureConst() const override
        {
            // Ignore the edge case in which all constituents are constant and
            // equal, as this is not common enough to warrant the increase in
            // efficiency gained by treating it as constant
            return false;
        }

        //- Harmonic mass fraction mixing
        declareScalarFieldFunctions(harmonicMixture);

        //- Return model dependence
        materialDepList();

        //- Update table
        virtual void updateTable(const word& modelName) override;

        //- Models casting
        virtual baseModels<scalar>* castScalarModel
        (
            const word& modelName
        ) override
        {
            if (harmonicMixtureModel::typeName == modelName)
            {
                return dynamic_cast<harmonicMixtureModel*>(this);
            }
            return nullptr;
        }

        //- Cast to tensor model
        virtual bool mixture() const override
        {
            return true;
        }

        //- Re-read data
        virtual bool read() override
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
