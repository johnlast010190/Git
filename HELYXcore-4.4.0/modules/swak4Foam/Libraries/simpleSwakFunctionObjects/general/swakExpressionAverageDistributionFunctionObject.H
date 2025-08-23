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
    (c) 2024 Engys Ltd.

Class
    swakExpressionAverageDistributionFunctionObject

Description
    Averages of a swak-expression in various directions

SourceFiles
    swakExpressionAverageDistributionFunctionObject.C

Contributors/Copyright:
    2008-2011, 2013-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef swakExpressionAverageDistributionFunctionObject_H
#define swakExpressionAverageDistributionFunctionObject_H

#include "timelineFunctionObject/timelineFunctionObject.H"

#include "CommonValueExpressionDriver.H"

#include "DataStructures/Distribution/SimpleDistribution.H"

#include "primitives/Vector/vector/vector.H"
#include "primitives/Tensor/tensor/tensor.H"
#include "primitives/SymmTensor/symmTensor/symmTensor.H"
#include "primitives/SphericalTensor/sphericalTensor/sphericalTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class swakExpressionAverageDistributionFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class swakExpressionAverageDistributionFunctionObject
:
    public timelineFunctionObject
{
    // Private Data

        //- The expression that is to be evaluated
        exprString expression_;

        //- The expression for the weight
        exprString weightExpression_;

        //- The weight values
        autoPtr<Field<scalar>> weightValues_;

        //- The expression for the mask
        exprString maskExpression_;

        //- The mask values
        autoPtr<Field<bool>> maskValues_;

        //- The expression for the abscissa
        exprString abscissaExpression_;

        //- The expression driver
        autoPtr<CommonValueExpressionDriver> driver_;

        //- Calculate extremes in the distribution from the current abscissa
        bool dynamicExtremesAbscissa_;


    // Private Member Functions

        template <typename T>
        void getDistributionInternal(PtrList<SimpleDistribution<T>>& dist);

        template <typename T>
        void writeADistribution(const PtrList<SimpleDistribution<T>>& dist);

        template <typename T>
        void reportADistribution(const PtrList<SimpleDistribution<T>>& dist);


protected:

    // Protected Data

        PtrList<SimpleDistribution<scalar>> distScalar_;
        PtrList<SimpleDistribution<vector>> distVector_;
        PtrList<SimpleDistribution<tensor>> distTensor_;
        PtrList<SimpleDistribution<symmTensor>> distSymmTensor_;
        PtrList<SimpleDistribution<sphericalTensor>> distSphericalTensor_;

        wordList componentNames_;


    // Protected Member Functions

        word dirName();

        virtual word baseName();

        template<class Type, class AType>
        autoPtr<PtrList<SimpleDistribution<AType>>> setData
        (
            const Field<AType>& xValues,
            const Field<Type>& values,
            const Field<scalar>& weights
        );

        template<class Type, class AType>
        autoPtr<PtrList<SimpleDistribution<AType>>> setData
        (
            const Field<AType>& xValues,
            const Field<Type>& values,
            const Field<scalar>& weights,
            const Field<bool>& mask
        );

        void clearDistributions();


public:

    //- Runtime type information
    TypeName("swakExpressionAverageDistribution");


    // Constructors

        //- Construct from components
        swakExpressionAverageDistributionFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        swakExpressionAverageDistributionFunctionObject
        (
            const swakExpressionAverageDistributionFunctionObject&
        ) = delete;


    // Member Functions

        virtual void writeSimple();

        virtual bool start();


protected:

    // Protected Member Functions

        virtual stringList columnNames();

        virtual wordList fileNames();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const swakExpressionAverageDistributionFunctionObject&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "swakExpressionAverageDistributionFunctionObjectTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
