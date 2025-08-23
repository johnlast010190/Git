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
    (c) 2016 Engys Ltd.
    (c) 2011-2023 OpenFOAM Foundation

Class
    Foam::icoPolynomial

Group
    grpSpecieEquationOfState

Description
    Incompressible, polynomial form of equation of state, using a polynomial
    function for density.

    Coefficient mixing is very inaccurate and not supported, so this equation of
    state is not applicable to mixtures.

    The polynomial expression is evaluated from:

    \verbatim
            rho = 1000 - 0.05*T + 0.003*sqr(T)
    \endverbatim

Usage
    \table
        Property     | Description
        rhoCoeffs    | Density polynomial coefficients
    \endtable

    Example specification of the icoPolynomial equation of state:
    \verbatim
    equationOfState
    {
        rhoCoeffs<8>    ( 1000 -0.05 0.003 0 0 0 0 0 );
    }
    \endverbatim

SourceFiles
    icoPolynomialI.H
    icoPolynomial.C

See also
    Foam::Polynomial

\*---------------------------------------------------------------------------*/

#ifndef icoPolynomial_H
#define icoPolynomial_H

#include "memory/autoPtr/autoPtr.H"
#include "primitives/functions/Polynomial/Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie, int PolySize>
class icoPolynomial;

template<class Specie, int PolySize>
icoPolynomial<Specie, PolySize> operator+
(
    const icoPolynomial<Specie, PolySize>&,
    const icoPolynomial<Specie, PolySize>&
);

template<class Specie, int PolySize>
icoPolynomial<Specie, PolySize> operator*
(
    const scalar,
    const icoPolynomial<Specie, PolySize>&
);

template<class Specie, int PolySize>
icoPolynomial<Specie, PolySize> operator==
(
    const icoPolynomial<Specie, PolySize>&,
    const icoPolynomial<Specie, PolySize>&
);

template<class Specie, int PolySize>
Ostream& operator<<
(
    Ostream&,
    const icoPolynomial<Specie, PolySize>&
);


/*---------------------------------------------------------------------------*\
                        Class icoPolynomial Declaration
\*---------------------------------------------------------------------------*/

template<class Specie, int PolySize=8>
class icoPolynomial
:
    public Specie
{
    // Private data

        //- Density polynomial coefficients [kg/m^3/K^i]
        Polynomial<PolySize> rhoCoeffs_;

        //- polynomial limiters
        scalar rhoMin_;
        scalar rhoMax_;


public:

    // Constructors

        //- Construct as copy
        icoPolynomial(const icoPolynomial&) = default;

        //- Construct from components
        inline icoPolynomial
        (
            const Specie& sp,
            const fvMesh& mesh,
            const Polynomial<PolySize>& rhoPoly,
            scalar rhoMin,
            scalar rhoMax
        );

        //- Construct from dictionary
        icoPolynomial(const objectRegistry& obr, const dictionary& dict);

        //- Construct as named copy
        inline icoPolynomial(const word& name, const icoPolynomial&);

        //- Construct and return a clone
        inline autoPtr<icoPolynomial> clone() const;

        // Selector from dictionary
        inline static autoPtr<icoPolynomial> New
        (
            const objectRegistry& obr,
            const dictionary& dict
        );


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "icoPolynomial<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = true;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy contribution [J/kg]
            inline scalar h(const scalar p, const scalar T) const;

            //- Return Cp contribution [J/(kg K]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return internal energy contribution [J/kg]
            inline scalar e(const scalar p, const scalar T) const;

            //- Return Cv contribution [J/(kg K]
            inline scalar Cv(scalar p, scalar T) const;

            //- Return entropy correction to the integral of Cp/T [J/kg/K]
            inline scalar sp(const scalar p, const scalar T) const;

            //- Return entropy correction to the integral of Cv/T [J/kg/K]
            inline scalar sv(const scalar p, const scalar T) const;

            //- Return compressibility [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            inline scalar CpMCv(scalar p, scalar T) const;

            //- Return volumetric coefficient of thermal expansion [1/T]
            inline scalar alphav(const scalar p, const scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators


        icoPolynomial& operator=(const icoPolynomial&) = default;
        inline void operator+=(const icoPolynomial&);
        inline void operator*=(const scalar);


    // Friend operators

        friend icoPolynomial operator+ <Specie, PolySize>
        (
            const icoPolynomial&,
            const icoPolynomial&
        );

        friend icoPolynomial operator* <Specie, PolySize>
        (
            const scalar s,
            const icoPolynomial&
        );

        friend icoPolynomial operator== <Specie, PolySize>
        (
            const icoPolynomial&,
            const icoPolynomial&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie, PolySize>
        (
            Ostream&,
            const icoPolynomial&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "equationOfState/icoPolynomial/icoPolynomialI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "equationOfState/icoPolynomial/icoPolynomial.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
