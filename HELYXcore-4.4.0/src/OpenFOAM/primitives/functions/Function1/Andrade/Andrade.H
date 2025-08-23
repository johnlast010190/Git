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
    (c) 2021 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

Class
    Foam::Function1Types::Andrade

Description
    Function one for transport properties using the Andrade function for
    the natural logarithm of dynamic viscosity
    and thermal conductivity of liquids:

    \verbatim
        log(mu) = coeffs[0] + coeffs[1]*T + coeffs[2]*sqr(T)
                + coeffs_[3]/(coeffs_[4] + T)
        log(kappa) = coeffs[0] + coeffs[1]*T + coeffs[2]*sqr(T)
                   + coeffs_[3]/(coeffs_[4] + T)
    );
    \endverbatim

    References:

    \verbatim
        Andrade, E. D. C. (1934).
        XLI. A theory of the viscosity of liquids.—Part I.
        The London, Edinburgh, and Dublin Philosophical Magazine
        and Journal of Science, 17(112), 497-511.
        Andrade, E. D. C. (1934).
        LVIII. A theory of the viscosity of liquids.—Part II.
        The London, Edinburgh, and Dublin Philosophical Magazine
        and Journal of Science, 17(113), 698-732.
    \endverbatim

Usage
    \table
        Property        | Description
        mu              | Dynamic viscosity polynomial coefficients
        kappa           | Thermal conductivity polynomial coefficients
    \endtable

    Example of the specification of the transport properties for water@200bar:
    \verbatim
    muModel scalarFunction1;
    muModelCoeffs
    {
        mu    (-25.8542 0.031256 -2.2e-05 3289.918 -11.4784);
    }
    kappaModel scalarFunction1;
    kappaModelCoeffs
    {
        kappa (-2.56543 0.008794 -9.8e-06 100.368 0);
    }
    \endverbatim

SourceFiles
    Andrade.C

\*---------------------------------------------------------------------------*/

#ifndef Andrade_H
#define Andrade_H

#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class Andrade Declaration
\*---------------------------------------------------------------------------*/

class Andrade
:
    public Function1<scalar>
{
    // Private data

        //- Equation coefficients
        FixedList<scalar, 5> coeffs_;


public:

    // Runtime type information
    TypeName("Andrade");


    // Constructors

        //- Construct from entry name and dictionary
        Andrade
        (
            const word& entryName,
            const dictionary& dict
        );

        //- Construct from components
        Andrade
        (
            const word& entryName,
            const FixedList<scalar, 5>&
        );

        //- Copy constructor
        Andrade(const Andrade& se);

        //- Construct and return a clone
        virtual tmp<Function1<scalar>> clone() const
        {
            return tmp<Function1<scalar>>(new Andrade(*this));
        }


    //- Destructor
    virtual ~Andrade();


    // Member Functions

        //- Return value for input variable
        scalar value(const scalar x) const;

        //- Derivative
        scalar derivative(const scalar x) const
        {
            NotImplemented;
        }

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Andrade&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
