/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : dev
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
    (c) ICE Stroemungsfoschungs GmbH
    (c) 1991-2008 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::AccumulationCalculation

Description
    Calculates the accumulations for a field. The purpose of this utilitiy
    class is to cache important intermediate results and thus make sure that
    - they are not calculated twice
    - parallel communication is minimized
    Also should it collect the implementation of all possible accumulations
    in one place

SourceFiles
    AccumulationCalculation.C

Contributors/Copyright:
    2013-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef AccumulationCalculation_H
#define AccumulationCalculation_H

#include "namedEnums/NumericAccumulationNamedEnum.H"
#include "CommonValueExpressionDriver.H"
#include "DataStructures/Distribution/SimpleDistribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class AccumulationCalculation Declaration
\*---------------------------------------------------------------------------*/

template <typename Type>
class AccumulationCalculation
{
    // Private data

        //- The data
        const Field<Type> &data_;

        //- Mask that selects what part of the data to use
        const Field<bool> mask_;

        //- The data that is actually used
        autoPtr<Field<Type>> usedData_;

        //- Is this point data?
        const bool isPoint_;

        //- The driver
        const CommonValueExpressionDriver &driver_;

        //- Weights on demand
        autoPtr<scalarField> weights_;

        //- On demand distribution
        autoPtr<SimpleDistribution<Type>> distribution_;

        //- On demand weighted distribution
        autoPtr<SimpleDistribution<Type>> weightedDistribution_;

        bool hasWeightSum_;
        scalar weightSum_;

        bool hasSize_;
        label size_;

        bool hasMaximum_;
        Type maximum_;

        bool hasMinimum_;
        Type minimum_;

        bool hasAverage_;
        Type average_;

        bool hasWeightedAverage_;
        Type weightedAverage_;

        bool hasSum_;
        Type sum_;

        bool hasWeightedSum_;
        Type weightedSum_;

        bool hasSumMag_;
        Type sumMag_;

        static const label defaultNumberOfBins_=100;

        label numberOfBins_;


    // Private Member Functions

        //- Construct null
        AccumulationCalculation(const IOobject&);

        label maskSize() const;


protected:

    const scalarField& weights();
    const Field<Type>& data();


public:

    // Constructors

        //- Constructor
        AccumulationCalculation
        (
            const Field<Type>& data,
            bool isPoint,
            CommonValueExpressionDriver& driver
        );

         //- Constructor with a mask
        AccumulationCalculation
        (
            const Field<Type>& data,
            bool isPoint,
            CommonValueExpressionDriver& driver,
            const Field<bool>& mask
        );

        //- Constructor with weights
        AccumulationCalculation
        (
            const Field<Type>& data,
            bool isPoint,
            CommonValueExpressionDriver& driver,
            const Field<scalar>& weight
        );

        //- Constructor with a mask and weights
        AccumulationCalculation
        (
            const Field<Type>& data,
            bool isPoint,
            CommonValueExpressionDriver& driver,
            const Field<bool>& mask,
            const Field<scalar>& weight
        );


    //- Destructor
    ~AccumulationCalculation();


    // Member Functions

        const SimpleDistribution<Type>& distribution();

        const SimpleDistribution<Type>& weightedDistribution();

        label size();
        scalar weightSum();

        Type maximum();
        Type minimum();
        Type average();
        Type weightedAverage();
        Type sum();
        Type sumMag();
        Type weightedSum();

        void resetNumberOfBins
        (
            const label newNumberOfBins,
            const scalar binWidth
        );


    // Member Operators

        //- Do an accumulation according to a specification
        Type operator()
        (
            const NumericAccumulationNamedEnum::accuSpecification&
        );

        //- Disallow default bitwise assignment
        void operator=(const AccumulationCalculation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "AccumulationCalculation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
