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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2019-2024 Engys Ltd.

Class
    Foam::Function1Types::TableBase

Description
    Base class for table with bounds handling, interpolation and integration

SourceFiles
    TableBase.C

\*---------------------------------------------------------------------------*/

#ifndef TableBase_H
#define TableBase_H

#include "primitives/functions/Function1/Function1/Function1.H"
#include "primitives/Tuple2/Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class interpolationWeights;

namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                        Class TableBase Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class TableBase
:
    public Function1<Type>
{
public:

    // Public data types

        //- Enumeration for handling out-of-bound values
        enum boundsHandling
        {
            ERROR,          //!< Exit with a FatalError
            WARN,           //!< Issue warning and clamp value (default)
            CLAMP,          //!< Clamp value to the start/end value
            REPEAT          //!< Treat as a repeating list
        };


protected:

    // Protected data

        //- Table name
        const word name_;

        //- Enumeration for handling out-of-bound values
        const boundsHandling boundsHandling_;

        //- Interpolation type
        const word interpolationScheme_;

        //- Table data
        mutable List<Tuple2<scalar, Type>> table_;

        //- Extracted values
        mutable autoPtr<scalarField> tableSamplesPtr_;

        //- Interpolator method
        mutable autoPtr<interpolationWeights> interpolatorPtr_;

        //- Cached indices and weights
        mutable labelList currentIndices_;

        mutable scalarField currentWeights_;


    // Protected Member Functions

        //- Return (demand driven) interpolator
        const interpolationWeights& interpolator() const;


public:

    // Constructors

        //- Construct from dictionary - note table is not populated
        TableBase(const word& name, const dictionary& dict);

        //- Copy constructor. Note: steals interpolator, tableSamples
        TableBase(const TableBase<Type>& tbl);


    //- Destructor
    virtual ~TableBase();


    // Member Functions

        //- Return the out-of-bounds handling as a word
        word boundsHandlingToWord(const boundsHandling& bound) const;

        //- Return the out-of-bounds handling as an enumeration
        boundsHandling wordToBoundsHandling(const word& bound) const;

        //- Set the out-of-bounds handling from enum, return previous setting
        boundsHandling outOfBounds(const boundsHandling& bound);

        //- Check the table for size and consistency
        void check() const;

        //- Check minimum table bounds
        bool checkMinBounds(const scalar x, scalar& xDash) const;

        //- Check maximum table bounds
        bool checkMaxBounds(const scalar x, scalar& xDash) const;

        //- Convert time
        virtual void convertTimeBase(const Time& t) override;

        //- Return Table value
        virtual Type value(const scalar x) const override;

        //- Integrate between two (scalar) values
        virtual Type integrate
        (
            const scalar x1,
            const scalar x2
        ) const override;

        //- Integrate function value divided by abscissa between two
        // (scalar) values
        virtual Type integrateYoverX
        (
            const scalar x1,
            const scalar x2
        ) const override;

        //- Return the derivative of the table
        virtual Type derivative(const scalar x) const override;

        //- Return the reference values
        virtual tmp<scalarField> x() const;

        //- Return the dependent values
        virtual tmp<Field<Type>> y() const;

        //- Write all table data in dictionary format
        virtual void writeData(Ostream& os) const override;

        //- Write keywords only in dictionary format. Used for non-inline
        //  table types
        virtual void writeEntries(Ostream& os) const;

        //- Shift time
        virtual void shiftTimeBase(const Time& t) override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const TableBase<Type>&) = delete;
};


// * * * * * * * * * Template specialisations  * * * * * * * * * * * * * * * //

//- Bool Type

template<>
bool TableBase<bool>::value(const scalar x) const;

template<>
bool TableBase<bool>::integrate(const scalar x1, const scalar x2) const;

template<>
bool TableBase<bool>::integrateYoverX(const scalar x1, const scalar x2) const;

template<>
bool TableBase<bool>::derivative(const scalar x) const;

template<>
Foam::scalarField TableBase<Foam::scalarField>::derivative(const scalar x) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "primitives/functions/Function1/Table/TableBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
