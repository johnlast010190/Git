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
    SimpleDistribution

Description
    Encapsulates the regular Distribution-class to allow fast usage with fields
    and parallelize

SourceFiles
    SimpleDistribution.C

Contributors/Copyright:
    2013-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef SimpleDistribution_H
#define SimpleDistribution_H

#include "Distribution.H"
#include "primitives/Tuple2/Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class SimpleDistribution;

template<class Type>
Istream& operator>>(Istream&, SimpleDistribution<Type>&);

template<class Type>
Ostream& operator<<(Ostream&, const SimpleDistribution<Type>&);

template<class Type>
SimpleDistribution<Type> operator+
(
    const SimpleDistribution<Type>&,
    const SimpleDistribution<Type>&
);

/*---------------------------------------------------------------------------*\
                    Class SimpleDistribution Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class SimpleDistribution
:
    protected Distribution<Type>
{
    // Private data

        //- Store the valid limits
        List<Pair<label>> validLimits_;

        bool hasInvalidValue_;

        Type invalidValue_;


   // Private Member Functions

        void recalcLimits();

        void setInvalidValue(const Type &val);


protected:

    // Protected Data

        //- Minimum for distributions that can handle it
        List<Field<scalar>> minimum_;

        //- Maximum for distributions that can handle it
        List<Field<scalar>> maximum_;

        //- Number of samples for distributions that can handle it
        List<Field<label>> nSamples_;


    // Protected Member Functions

        //- Returns the indices of the first and last non-zero entries
        Pair<label> validLimits(direction cmpt) const;

        //- round value down to the next nearst power of 10
        Type smaller10Power(const Type value);


public:

    // Constructors

        //- Construct null
        SimpleDistribution();

        //- Setting a size
        SimpleDistribution(const scalar binSize);

        //- Setting using a span plus the bin number
        SimpleDistribution(const Type min, const Type max, const label binNr);

        //- From a base class
        SimpleDistribution(const Distribution<Type>&);

        //- Direct copy
        SimpleDistribution(const SimpleDistribution<Type>&);


    // Member Functions

        void calc
        (
            const Field<Type>& values,
            const Field<Type>& weights
        );

        void calc
        (
            const Field<Type>& values,
            const Field<Type>& weights,
            const Field<bool>& mask
        );

       void calcScalarWeight
       (
            const Field<Type>& values,
            const Field<scalar>& weights,
            bool doReduce = true
        );

       void calcScalarWeight
       (
            const Field<Type>& values,
            const Field<scalar>& weights,
            const Field<bool>& mask,
            bool doReduce = true
        );

       void calcMinimumMaximum
       (
            const Field<Type>& values,
            const Field<scalar>& weights,
            const Field<bool>& mask
        );

       void calcMinimumMaximum
       (
            const Field<Type>& values,
            const Field<scalar>& weights
        );

        void divideByDistribution
        (
            const SimpleDistribution<Type>& weightSum,
            const Type& valueIfZero
        );

        //- Return the distribution of the total bin weights - adapted to allow negative weights
        List<List<Pair<scalar>>> rawNegative() const;

        template<class FType>
        List<List<Tuple2<scalar, FType>>> rawField
        (
            const List<Field<FType>>& f
        ) const;

        Type mean() const
        {
            return Distribution<Type>::mean();
        }

        Type median() const
        {
            return Distribution<Type>::median();
        }

        scalar totalWeight(direction cmpt) const
        {
            return Distribution<Type>::totalWeight(cmpt);
        }

        Type min() const;

        Type max() const;

        //- Return the value where a fraction of the distribution is smaller
        Type quantile(scalar frac) const;

        //- Return the fraction of the distribution that is smaller than the value
        Type smaller(scalar value) const;

        void write(const fileName& filePrefix) const
        {
            Distribution<Type>::write(filePrefix);
        }

        label maxNrBins() const;

        //- Write the distribution to file: only key and raw.
        //  Produces a separate file for each component.
        void writeRaw(const fileName& filePrefix) const;


    // Member Operators

        void operator=(const SimpleDistribution<Type>&);

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const SimpleDistribution<Type>&
        );

        friend Istream& operator>> <Type>
        (
            Istream& is,
            SimpleDistribution<Type>& d
        );

        friend SimpleDistribution<Type> operator+ <Type>
        (
            const SimpleDistribution<Type>& d1,
            const SimpleDistribution<Type>& d2
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#   include "SimpleDistribution.C"
#endif


#endif

// ************************************************************************* //
