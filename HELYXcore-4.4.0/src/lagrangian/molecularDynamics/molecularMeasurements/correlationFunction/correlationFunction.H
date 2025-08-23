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
    (c) 2024 Engys Ltd.

Class
    Foam::correlationFunction

Description

SourceFiles
    correlationFunctionI.H
    correlationFunction.C
    correlationFunctionIO.C

\*---------------------------------------------------------------------------*/

#ifndef correlationFunction_H
#define correlationFunction_H

#include "bufferedAccumulator/bufferedAccumulator.H"
#include "db/dictionary/dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class correlationFunction;

template<class Type>
Ostream& operator<<
(
    Ostream&,
    const correlationFunction<Type>&
);

/*---------------------------------------------------------------------------*\
                    Class correlationFunction Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class correlationFunction
:
    public bufferedAccumulator<scalar>
{
    // Private data

        const polyMesh& mesh_;

        Field<Field<Type>> tZeroBuffers_;

        scalar duration_;
        scalar sampleInterval_;
        scalar averagingInterval_;

        label sampleSteps_;


    // Private Member Functions

        void setTimesAndSizes(const label);


public:

    //- Component type
    typedef typename pTraits<Type>::cmptType cmptType;


    // Static data members

        static const char* const typeName;


    // Constructors

        //- Construct from dictionary
        correlationFunction
        (
            const polyMesh&,
            const dictionary&,
            const label tZeroBufferSize
        );

        //- Construct from components
        correlationFunction
        (
            const polyMesh&,
            const label tZeroBufferSize,
            const scalar duration,
            const scalar sampleInterval,
            const scalar averagingInterval

        );

        //- Disallow default bitwise copy construction
        correlationFunction(const correlationFunction<Type>&) = delete;


    //- Destructor
    ~correlationFunction();


    // Member Functions

        void calculateCorrelationFunction(const Field<Type>&);

        void calculateCorrelationFunction(const Type&);

        scalar integral() const;

        bool writeAveraged(Ostream&) const;


    // Access

        inline const Field<Field<Type>>& tZeroBuffers() const;

        inline scalar duration() const;

        inline scalar sampleInterval() const;

        inline scalar averagingInterval() const;

        inline label sampleSteps() const;

        inline label measurandFieldSize() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const correlationFunction<Type>&) = delete;


    // IOstream Operators

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const correlationFunction<Type>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "correlationFunction/correlationFunctionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "correlationFunction/correlationFunction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
