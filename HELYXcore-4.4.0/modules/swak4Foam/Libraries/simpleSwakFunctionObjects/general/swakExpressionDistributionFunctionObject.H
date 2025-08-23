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

Class
    swakExpressionDistributionFunctionObject

Description
    Distribution of a swak-expression

SourceFiles
    swakExpressionDistributionFunctionObject.C

Contributors/Copyright:
    2008-2011, 2013-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef swakExpressionDistributionFunctionObject_H
#define swakExpressionDistributionFunctionObject_H

#include "distributionFunctionObject/distributionFunctionObject.H"

#include "CommonValueExpressionDriver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class swakExpressionDistributionFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class swakExpressionDistributionFunctionObject
:
    public distributionFunctionObject
{
    // Private Data

        //- The expression that is to be evaluated
        exprString expression_;

        //- The expression for the weight
        exprString weightExpression_;

        //- The expression for the mask
        exprString maskExpression_;

        //- The expression driver
        autoPtr<CommonValueExpressionDriver> driver_;

        autoPtr<Field<scalar>> weightValuesScalar_;
        autoPtr<Field<vector>> weightValuesVector_;
        autoPtr<Field<tensor>> weightValuesTensor_;
        autoPtr<Field<symmTensor>> weightValuesSymmTensor_;
        autoPtr<Field<sphericalTensor>> weightValuesSphericalTensor_;
        autoPtr<Field<bool>> maskValues_;


    // Private Member Functions

        template <typename T>
        void getDistributionInternal
        (
            autoPtr<SimpleDistribution<T>>& dist,
            autoPtr<Field<T>>& sameWeight
        );


protected:

    // Protected Member Functions

        word dirName();

        virtual word baseName();

        virtual void getDistribution();


public:

    //- Runtime type information
    TypeName("swakExpressionDistribution");


    // Constructors

        //- Construct from components
        swakExpressionDistributionFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        swakExpressionDistributionFunctionObject
        (
            const swakExpressionDistributionFunctionObject&
        ) = delete;


    // Member Functions

        void writeSimple();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const swakExpressionDistributionFunctionObject&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "swakExpressionDistributionFunctionObjectTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
