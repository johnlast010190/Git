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
    (c) held by original author

Class
    BisectionRoot

Description
    Bisection root Based on Numerical Recipes in C++, Section 9.1, page 358.
    Function is provided as a template parameter function object, evaluated
    using operator()(const scalar x)

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    BisectionRoot.C

\*---------------------------------------------------------------------------*/

#ifndef BisectionRoot_H
#define BisectionRoot_H

#include "primitives/Scalar/scalar/scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class BisectionRoot Declaration
\*---------------------------------------------------------------------------*/

template<class Func>
class BisectionRoot
{
private:

    // Private data

        //- Reference to a function
        const Func& f_;

        //- Required accuracy
        const scalar eps_;


public:

    // Static data members

        //- Maximum number of iterations
        static const label maxIter;


    // Constructors

        //- Construct given a function
        BisectionRoot(const Func& f, const scalar eps);


    // Destructor - default


    // Member Functions

        //- Return root
        scalar root(const scalar x0, const scalar x1) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "findroot/BisectionRoot.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
