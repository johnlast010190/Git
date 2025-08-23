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
    distributionFunctionObject

Description
    Basis for function objects that generate distributions

SourceFiles
    distributionFunctionObject.C

Contributors/Copyright:
    2008-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef distributionFunctionObject_H
#define distributionFunctionObject_H

#include "timelineFunctionObject/timelineFunctionObject.H"
#include "primitives/enums/NamedEnum.H"

#include "DataStructures/Distribution/SimpleDistribution.H"

#include "primitives/Vector/vector/vector.H"
#include "primitives/Tensor/tensor/tensor.H"
#include "primitives/SymmTensor/symmTensor/symmTensor.H"
#include "primitives/SphericalTensor/sphericalTensor/sphericalTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class distributionFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class distributionFunctionObject
:
    public timelineFunctionObject
{
    // Private Data

        bool writeTimeline_;

        bool writeDistribution_;

        bool startup_;
        word startTime_;


    // Private Member Functions

        template<class Type>
        void writeATimeline(const autoPtr<SimpleDistribution<Type>>& pDist);

        template<class Type>
        void writeADistribution
        (
            const autoPtr<SimpleDistribution<Type>>& pDist
        );

        template<class Type>
        void reportADistribution
        (
            const autoPtr<SimpleDistribution<Type>>& pDist
        );


protected:

    // Protected Data

        autoPtr<SimpleDistribution<scalar>> distScalar_;
        autoPtr<SimpleDistribution<vector>> distVector_;
        autoPtr<SimpleDistribution<tensor>> distTensor_;
        autoPtr<SimpleDistribution<symmTensor>> distSymmTensor_;
        autoPtr<SimpleDistribution<sphericalTensor>> distSphericalTensor_;


    // Protected Member Functions

        template<class Type>
        autoPtr<SimpleDistribution<Type>> setData
        (
            const Field<Type>& values,
            const Field<Type>& weights,
            const scalar binWidth = -1
        );

        template<class Type>
        autoPtr<SimpleDistribution<Type>> setData
        (
            const Field<Type>& values,
            const Field<Type>& weights,
            const Field<bool>& mask,
            const scalar binWidth = -1
        );

        template<class Type>
        autoPtr<SimpleDistribution<Type>> setDataScalar
        (
            const Field<Type>& values,
            const Field<scalar>& weights,
            const scalar binWidth = -1
        );

        template<class Type>
        autoPtr<SimpleDistribution<Type>> setDataScalar
        (
            const Field<Type>& values,
            const Field<scalar>& weights,
            const Field<bool>& mask,
            const scalar binWidth = -1
        );

        void clearDistributions();

        virtual fileName dataDir();


public:

    //- Runtime type information
    TypeName("distributionFunctionObject");


    // Constructors

        //- Construct from components
        distributionFunctionObject
        (
            const word &name,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        distributionFunctionObject(const distributionFunctionObject&) = delete;


    // Member Functions

        bool writeTimeline() const
        {
            return writeTimeline_;
        }

        bool writeDistribution() const
        {
            return writeDistribution_;
        }

        virtual void flush();

        virtual void writeSimple();

        virtual bool start();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const distributionFunctionObject&) = delete;


protected:

    // Protected Member Functions

        virtual stringList columnNames();

        virtual wordList fileNames();

        virtual word baseName() = 0;

        virtual stringList componentNames();

        virtual void getDistribution() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "distributionFunctionObjectI.H"
#endif

#endif

// ************************************************************************* //
