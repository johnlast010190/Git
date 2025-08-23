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
    (c) 2011-2024 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

Class
    Foam::Function1Types::NonUniformTable

Description
    Non-uniform tabulated property function that linearly interpolates between
    the values.

    To speed-up the search of the non-uniform table a uniform jump-table is
    created on construction which is used for fast indirect addressing into
    the table.

Usage
    \table
        Property    | Description
        values      | List of value pairs
    \endtable

    Example for the density of water between 280 and 350K
    \verbatim
    rho
    {
        type    nonUniformTable;

        values
        (
            (280 999.87)
            (300 995.1)
            (350 973.7)
        );
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef NonUniformTable1_H
#define NonUniformTable1_H

#include "primitives/functions/Function1/Function1/Function1.H"
#include "primitives/Tuple2/Tuple2.H"
#include "interpolations/interpolationTable/tableReaders/tableReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                       Class NonUniformTable Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class NonUniformTable
:
    public Function1<Type>
{
    // Private member data

        //- Lowest value in the NonUniformTable
        scalar low_;

        //- Highest value in the NonUniformTable
        scalar high_;

        //- Table values
        List<Tuple2<scalar, Type>> values_;

        //- Increment derived from low_, high_ and values_.size()
        scalar delta_;

        //- The jump table
        List<label> jumpTable_;


protected:

        //- Return the lower index of the interval in the table
        //  corresponding to the given value
        inline label index(scalar x) const;


public:

    //- Runtime type information
    TypeName("nonUniformTable");


    // Constructors

        //- Construct from name and dictionary
        NonUniformTable(const word& name, const dictionary& dict);

        //- Copy constructor
        NonUniformTable(const NonUniformTable<Type>& nut);


    // Member Functions

        //- Return the non-uniform table of values
        const List<Tuple2<scalar, Type>>& values() const
        {
            return values_;
        }

        //- Evaluate the function and return the result
        virtual Type value(scalar x) const;

        //- Integrate between two scalar values
        scalar integrate(const scalar x) const;

        //- Integrate the function divided by x and return the result
        scalar integrateByX(const scalar x) const;

        //- Integrate the function divided by x and return the result
        scalar integrateByX(const scalar x1, const scalar x2) const;

        //- Integrate between two scalar values
        Type integrate(const scalar x1, const scalar x2) const;

        //- Evaluate the derivative of the function and return the result
        Type derivative(scalar T) const;

        //- Write the function coefficients
        void writeData(Ostream& os) const;


    // Member Operators

        //- Assignment
        void operator=(const NonUniformTable<Type>&);
};


template<>
scalar NonUniformTable<scalar>::integrate
(
    const scalar x1,
    const scalar x2
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "primitives/functions/Function1/NonUniformTable/NonUniformTable1I.H"

#ifdef NoRepository
    #include "primitives/functions/Function1/NonUniformTable/NonUniformTable1.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
