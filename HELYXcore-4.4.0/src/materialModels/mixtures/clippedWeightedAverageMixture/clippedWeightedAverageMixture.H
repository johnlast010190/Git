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
    Foam::clippedWeightedAverageMixture

Description
    Homogeneus mixture model.

SourceFiles
    clippedWeightedAverageMixture.C

\*---------------------------------------------------------------------------*/

#ifndef clippedWeightedAverageMixture_H
#define clippedWeightedAverageMixture_H

#include "memory/autoPtr/autoPtr.H"
#include "materialModels/materialModel/materialModel.H"
#include "mixtures/speciesMassFractions/speciesMassFractions.H"
#include "mixtures/phaseVolumeFractions/phaseVolumeFractions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class clippedWeightedAverageMixture Declaration
\*---------------------------------------------------------------------------*/

class clippedWeightedAverageMixture
:
    public materialModel,
    public clippedWeightedAverageMixtureModel,
    public vectorClippedWeightedAverageMixtureModel,
    public tensorClippedWeightedAverageMixtureModel
{
    // Private member data

        //- Reference to the class providing mixing fractions
        const volumeMassFractions& frac_;

        //- Lookup or construct fractions class
        const volumeMassFractions& lookupOrConstructBase();

public:

    // Constructors

        //- Construct from dictionary
        clippedWeightedAverageMixture
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& name,
            const word& phaseName,
            const word& specieName
        );

        //- New constructor
        static autoPtr<clippedWeightedAverageMixture> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );



    //- Destructor
    virtual ~clippedWeightedAverageMixture()
    {}

    // Member functions

        //- Return the instantiated type name
        TypeName("clippedWeightedAverageMixture");


    // Member functions

        //- Internal field mixing
        template<class Type>
        tmp<Field<Type>> clippedWeightedAverageMixtureTypeInternal
        (
            const UPtrList<baseModels<Type>>& mod
        ) const;

        //- Geometric field mixing
        template<class Type>
        tmp<VolField<Type>>
        clippedWeightedAverageMixtureTypeGeometric
        (
            const UPtrList<baseModels<Type>>& mod
        ) const;

        //- Patch mixing
        template<class Type>
        tmp<Field<Type>> clippedWeightedAverageMixtureTypePatch
        (
            const UPtrList<baseModels<Type>>& mod,
            const label patchi
        ) const;

        //- Cell mixing
        template<class Type>
        Type clippedWeightedAverageMixtureTypeCell
        (
            const UPtrList<baseModels<Type>>& mod,
            const label celli
        ) const;

        bool isclippedWeightedAverageMixtureConst() const override
        {
            // Ignore the edge case in which all constituents are constant and
            // equal, as this is not common enough to warrant the increase in
            // efficiency gained by treating it as constant
            return false;
        }

        bool isvectorClippedWeightedAverageMixtureConst() const override
        {
            return false;
        }

        bool istensorClippedWeightedAverageMixtureConst() const override
        {
            return false;
        }

        declareScalarFieldFunctions(clippedWeightedAverageMixture);

        declareVectorFieldFunctions(vectorClippedWeightedAverageMixture);

        declareTensorFieldFunctions(tensorClippedWeightedAverageMixture);

        //- Return model dependence
        materialDepList();

        //- Update table
        virtual void updateTable(const word& modelName) override;

        //- Cast to scalar model
        virtual baseModels<scalar>* castScalarModel
        (
            const word& modelName
        ) override;

        //- Cast to vector model
        virtual baseModels<vector>* castVectorModel
        (
            const word& modelName
        ) override;

        //- Cast to tensor model
        virtual baseModels<tensor>* castTensorModel
        (
            const word& modelName
        ) override;

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

#ifdef NoRepository
    #include "clippedWeightedAverageMixtureTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
