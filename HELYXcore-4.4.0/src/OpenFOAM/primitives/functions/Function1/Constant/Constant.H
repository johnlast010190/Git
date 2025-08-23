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
    (c) 2015 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2019-2022 Engys Ltd.

Class
    Foam::Function1Types::Constant

Description
    Templated function that returns a constant value.

    Usage - for entry \<entryName\> returning the value <value>:
    \verbatim
        <entryName>    constant  <value>
    \endverbatim

SourceFiles
    Constant.C

\*---------------------------------------------------------------------------*/

#ifndef Constant_H
#define Constant_H

#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class Constant Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Constant
:
    public Function1<Type>
{
    // Private data

        //- Constant value
        Type value_;


public:

    // Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from components
        Constant(const word& entryName, const Type& value);

        //- Construct from entry name and dictionary
        Constant(const word& entryName, const dictionary& dict);

        //- Construct from entry name and Istream
        //  Reads the constant value without the Function1 type
        //  for backward compatibility
        Constant(const word& entryName, Istream& is);

        //- Copy constructor
        Constant(const Constant<Type>& cnst);

        //- Construct and return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return tmp<Function1<Type>>(new Constant<Type>(*this));
        }


    //- Destructor
    virtual ~Constant();


    // Member Functions

        //- Return constant value
        Type value(const scalar) const;

        //- Integrate between two values
        Type integrate(const scalar x1, const scalar x2) const;

        //- Integrate function value divided by abscissa between two values
        Type integrateYoverX(const scalar x1, const scalar x2) const;

        //- Derivative
        Type derivative(const scalar x) const;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Constant<Type>&) = delete;
};


// * * * * * * * * * Template specialisations  * * * * * * * * * * * * * * * //

template<>
bool Constant<bool>::integrate(const scalar x1, const scalar x2) const;

template<>
bool Constant<bool>::integrateYoverX(const scalar x1, const scalar x2) const;

template<>
Foam::scalarField Foam::Function1Types::Constant<Foam::scalarField>::derivative
(
    const scalar x
) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "primitives/functions/Function1/Constant/Constant.C"
    #include "primitives/functions/Function1/Function1/Function1New.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
